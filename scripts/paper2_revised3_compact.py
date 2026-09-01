from pathlib import Path
import argparse
import zipfile,time,warnings,itertools
import numpy as np,pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score,f1_score
warnings.filterwarnings('ignore')
ROOT=Path('/mnt/data'); ZIP=ROOT/'ZENODO_UPLOAD_v1.0.zip'; OUT=ROOT/'paper2_revised3'; OUT.mkdir(exist_ok=True)
FEATURES=['voltage_V','current_A','power_W','frequency_Hz','power_factor','PE_m3','PE_m4','PE_m5','WPE','MSPE','CV','delta_PF','SII']; IDX={f:i for i,f in enumerate(FEATURES)}
FILES={'S1':'station1_normal_features.csv','S2':'station2_medium_features.csv','S3':'station3_high_features.csv'}; SPARSE={'inrush_surge','spike_transient','overload'}
SC={'VOLTAGE_RAMP':('mult',['voltage_V'],.10),'LOAD_RAMP':('mult',['current_A','power_W'],.20),'SENSOR_GAIN':('mult',['voltage_V','current_A','power_W'],.10),'FREQUENCY_OFFSET':('add',['frequency_Hz'],.50)}
BLOCK=3600; PERSIST=3; SEED=42; MIN_KNOWN_ACCEPT=.50

def stats(a):
 c=np.nanmedian(a,0); q1=np.nanquantile(a,.25,axis=0); q3=np.nanquantile(a,.75,axis=0); s=q3-q1; mad=1.4826*np.nanmedian(np.abs(a-c),0); sd=np.nanstd(a,0); s=np.where((np.isfinite(s))&(s>1e-12),s,mad); s=np.where((np.isfinite(s))&(s>1e-12),s,sd); s=np.where((np.isfinite(s))&(s>1e-12),s,1.0); return c,s

def load_raw(st):
 cols=FEATURES+['segment_id','anomaly','anomaly_type','event_id']
 with zipfile.ZipFile(ZIP) as z:
  with z.open('controlled-fault-benchmark-v1.0/'+FILES[st]) as f:return pd.read_csv(f,usecols=cols)

def load_static(): return {s:pd.read_csv(ROOT/'paper2_step1_adaptive_rlr'/f'static_{s}.csv') for s in FILES}

def thresholds(refX):
 c,s=stats(refX); ds=[]
 for q in range(0,len(refX),BLOCK):
  xb=refX[q:q+BLOCK]
  if len(xb)>=100: fc,_=stats(xb); ds.append(float(np.nanmax(np.abs((fc-c)/s))))
 ds=np.array(ds); gon=max(float(np.quantile(ds,.95)),.5); goff=max(float(np.quantile(ds,.75)),.70*gon); goff=.85*gon if goff>=gon else goff; return c,s,gon,goff

def drift(X,sc):
 kind,fs,amp=SC[sc]; Y=X.copy(); frac=np.linspace(0,1,len(X))
 for f in fs:
  j=IDX[f]; Y[:,j]=X[:,j]*(1+amp*frac) if kind=='mult' else X[:,j]+amp*frac
 return Y

def precompute(df,sc,zvals=(3.,4.,5.)):
 first=sorted(df.segment_id.dropna().unique())[0]; ref=df[df.segment_id.eq(first)]; ev=df[~df.segment_id.eq(first)].reset_index(drop=True)
 c0,s0,gon,goff=thresholds(ref[FEATURES].to_numpy(float)); X=drift(ev[FEATURES].to_numpy(float),sc); an=ev.anomaly.to_numpy(int); ids=ev.event_id.to_numpy(); ty=ev.anomaly_type.astype(object).to_numpy()
 blocks=[]
 for start in range(0,len(X),BLOCK):
  end=min(start+BLOCK,len(X)); xb=X[start:end]; fc,fs=stats(xb); cand={}
  dist=np.nanmax(np.abs((xb-fc)/fs),axis=1)
  for z in zvals:
   keep=np.isfinite(dist)&(dist<=z)
   if keep.sum()>=100: bc,bs=stats(xb[keep])
   else: bc=bs=None
   cand[z]=(bc,bs,int(keep.sum()),int(an[start:end][keep].sum()))
  m=an[start:end]==1
  blocks.append((xb,fc,fs,cand,ids[start:end][m],ty[start:end][m],xb[m]))
 return c0,s0,gon,goff,blocks

def simulate(pre,zgate=4.,eta=.20,mode='DT'):
 c,s,gon,goff,blocks=pre; c=c.copy(); s=s.copy(); rows=[]; active=(mode=='DUAL_ALWAYS'); onc=offc=0; acc=aa=act=trig=0
 t0=time.perf_counter()
 for xb,fc,fs,cand,eids,etypes,xanom in blocks:
  if len(eids):
   zz=(xanom-c)/s; a=pd.DataFrame(zz,columns=FEATURES); a.insert(0,'anomaly_type',etypes); a.insert(0,'event_id',eids); rows.append(a)
  D=float(np.nanmax(np.abs((fc-c)/s)))
  was_active=active
  if mode in ('DT','TRIGGER_ONLY') and not active:
   onc=onc+1 if D>gon else 0
   if onc>=PERSIST: active=True;trig+=1;onc=offc=0
  do_update = active if mode=='DUAL_ALWAYS' else was_active
  if do_update:
   act+=1
   if mode=='TRIGGER_ONLY': bc,bs=fc,fs; nk=len(xb); nab=int(np.sum([len(eids)]))
   else: bc,bs,nk,nab=cand[zgate]
   if bc is not None: c=(1-eta)*c+eta*bc; s=(1-eta)*s+eta*bs
   acc+=nk; aa+=nab
   if mode in ('DT','TRIGGER_ONLY'):
    Dpost=float(np.nanmax(np.abs((fc-c)/s)));offc=offc+1 if Dpost<goff else 0
    if offc>=PERSIST:active=False;onc=offc=0
 rep=pd.concat(rows,ignore_index=True).groupby(['event_id','anomaly_type'],as_index=False)[FEATURES].median()
 return rep,{'gamma_on':gon,'gamma_off':goff,'contamination':aa/acc if acc else np.nan,'accepted_update_rows':acc,'active_blocks':act,'trigger_count':trig,'runtime_seconds_simulation':time.perf_counter()-t0}

def impute(tr,te): med=tr[FEATURES].median();a=tr.copy();b=te.copy();a[FEATURES]=a[FEATURES].fillna(med);b[FEATURES]=b[FEATURES].fillna(med);return a,b

def fit(tr,n=20): m=RandomForestClassifier(n_estimators=n,class_weight='balanced_subsample',random_state=SEED,n_jobs=1,max_features='sqrt');m.fit(tr[FEATURES],tr.anomaly_type);return m

def ep(m,te): P=m.predict_proba(te[FEATURES]);e=-(P*np.log(np.clip(P,1e-12,1))).sum(1)/np.log(P.shape[1]);cls=np.asarray(m.classes_,object);return e,cls[P.argmax(1)]
def metrics(y,pred,score,tau,hs,known):
 y=np.asarray(y,object);u=np.isin(y,list(hs));k=~u;rej=score>tau;out=pred.copy();out[rej]='UNKNOWN';tar=y.copy();tar[u]='UNKNOWN';ka=np.mean(~rej[k]);kc=np.mean(out[k]==tar[k]);ur=np.mean(rej[u]);h=0 if kc+ur==0 else 2*kc*ur/(kc+ur);auc=roc_auc_score(u.astype(int),score);mf=f1_score(tar,out,labels=list(known)+['UNKNOWN'],average='macro',zero_division=0);return ka,kc,ur,h,auc,mf

def meta_tau(static,src,hs):
 hs=set(hs); pseudos=sorted(set(pd.concat([static[s] for s in src]).anomaly_type)-hs); rr=[];a,b=src
 for s,t in [(a,b),(b,a)]:
  for pseudo in pseudos:
   tr=static[s][~static[s].anomaly_type.isin(hs|{pseudo})].copy(); te=static[t][~static[t].anomaly_type.isin(hs)].copy();tr,te=impute(tr,te);m=fit(tr,10);score,_=ep(m,te);rr.append(pd.DataFrame({'score':score,'y':(te.anomaly_type.to_numpy()==pseudo).astype(int)}))
 q=pd.concat(rr);score=q.score.to_numpy();y=q.y.to_numpy();best=None
 for tau in np.unique(np.quantile(score,np.linspace(.40,.995,120))):
  ka=np.mean(score[y==0]<=tau)
  if ka<MIN_KNOWN_ACCEPT:continue
  ur=np.mean(score[y==1]>tau);h=0 if ka+ur==0 else 2*ka*ur/(ka+ur);key=(ur,h,-tau)
  if best is None or key>best[0]:best=(key,float(tau))
 return best[1]

static=load_static();main13=sorted(set().union(*[set(d.anomaly_type) for d in static.values()])-SPARSE)
for s in static:static[s]=static[s][static[s].anomaly_type.isin(main13)].reset_index(drop=True)
taus=pd.read_csv(ROOT/'paper2_step1_adaptive_rlr/step1_all_scenarios.csv').groupby(['test_station','heldout_type'])['tau'].first().to_dict()
ap=argparse.ArgumentParser(); ap.add_argument('--station',required=True); ap.add_argument('--mode',choices=['sens','multi','all'],default='all'); args=ap.parse_args();
all_sens=[];aud=[];multi=[];ab=[]
for tgt in [args.station]:
 print('station',tgt,flush=True);df=load_raw(tgt);src=[s for s in FILES if s!=tgt]
 models={}
 for held in main13:
  tr=pd.concat([static[s] for s in src]);tr=tr[tr.anomaly_type!=held].copy();tr,_=impute(tr,tr);models[held]=(tr,fit(tr),sorted(tr.anomaly_type.unique()),float(taus[(tgt,held)]))
 pre={sc:precompute(df,sc) for sc in SC}
 base={}
 for sc in SC:
  for z in [3.,4.,5.]:
   rep,a=simulate(pre[sc],z,.20,'DT'); rep=rep[rep.anomaly_type.isin(main13)].reset_index(drop=True); a.update({'test_station':tgt,'experiment':'zgate','setting':z,'scenario':sc});aud.append(a)
   if z==4:base[sc]=rep
   for held,(tr,m,known,tau) in models.items():
    _,te=impute(tr,rep);score,pred=ep(m,te);v=metrics(te.anomaly_type,pred,score,tau,{held},known);all_sens.append({'test_station':tgt,'heldout_type':held,'experiment':'zgate','setting':z,'scenario':sc,'known_acceptance':v[0],'known_correct_rate':v[1],'unknown_recall':v[2],'h_score':v[3],'unknown_auroc':v[4],'open_macro_f1':v[5]})
  for eta in [.10,.20,.40]:
   rep,a=simulate(pre[sc],4.,eta,'DT'); rep=rep[rep.anomaly_type.isin(main13)].reset_index(drop=True); a.update({'test_station':tgt,'experiment':'eta','setting':eta,'scenario':sc});aud.append(a)
   for held,(tr,m,known,tau) in models.items():
    _,te=impute(tr,rep);score,pred=ep(m,te);v=metrics(te.anomaly_type,pred,score,tau,{held},known);all_sens.append({'test_station':tgt,'heldout_type':held,'experiment':'eta','setting':eta,'scenario':sc,'known_acceptance':v[0],'known_correct_rate':v[1],'unknown_recall':v[2],'h_score':v[3],'unknown_auroc':v[4],'open_macro_f1':v[5]})
  rep,a=simulate(pre[sc],4.,.20,'TRIGGER_ONLY'); rep=rep[rep.anomaly_type.isin(main13)].reset_index(drop=True)
  for held,(tr,m,known,tau) in models.items():
   _,te=impute(tr,rep);score,pred=ep(m,te);v=metrics(te.anomaly_type,pred,score,tau,{held},known);ab.append({'test_station':tgt,'heldout_type':held,'scenario':sc,'method':'TRIGGER_ONLY_NO_GATE','known_acceptance':v[0],'known_correct_rate':v[1],'unknown_recall':v[2],'h_score':v[3],'unknown_auroc':v[4],'open_macro_f1':v[5]})
 print('sens/ab done',tgt,flush=True)
 pd.DataFrame(all_sens).to_csv(OUT/f'sensitivity_{tgt}.csv',index=False); pd.DataFrame(aud).to_csv(OUT/f'sensitivity_audit_{tgt}.csv',index=False); pd.DataFrame(ab).to_csv(OUT/f'ablation_trigger_only_{tgt}.csv',index=False)
 if args.mode=='sens': continue
 n=len(main13);starts=[0,3,6,9]; groups={2:[tuple(sorted([main13[i],main13[(i+1)%n]])) for i in starts],3:[tuple(sorted([main13[i],main13[(i+1)%n],main13[(i+2)%n]])) for i in starts]}
 for k,gg in groups.items():
  for hs in gg:
   hs=set(hs);tr=pd.concat([static[s] for s in src]);tr=tr[~tr.anomaly_type.isin(hs)].copy();tr,_=impute(tr,tr);m=fit(tr);known=sorted(tr.anomaly_type.unique());tau=meta_tau(static,src,hs)
   for sc in ['LOAD_RAMP','FREQUENCY_OFFSET']:
    _,te=impute(tr,base[sc]);score,pred=ep(m,te);v=metrics(te.anomaly_type,pred,score,tau,hs,known);multi.append({'test_station':tgt,'n_unknown_classes':k,'heldout_types':'|'.join(sorted(hs)),'scenario':sc,'known_acceptance':v[0],'known_correct_rate':v[1],'unknown_recall':v[2],'h_score':v[3],'unknown_auroc':v[4],'open_macro_f1':v[5],'tau':tau})
 print('multi done',tgt,flush=True)
 pd.DataFrame(all_sens).to_csv(OUT/f'sensitivity_{tgt}.csv',index=False); pd.DataFrame(aud).to_csv(OUT/f'sensitivity_audit_{tgt}.csv',index=False); pd.DataFrame(multi).to_csv(OUT/f'multi_unknown_{tgt}.csv',index=False); pd.DataFrame(ab).to_csv(OUT/f'ablation_trigger_only_{tgt}.csv',index=False)
print('done all')
