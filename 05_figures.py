"""Stage 5: poster figures from the JSON outputs (300 dpi PNG)."""
import json, glob, os, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
NAVY="#03163C"; TEAL="#0E7C86"; GOLD="#D99A00"; GREY="#8A94A6"; INK="#1B2333"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":15,"axes.edgecolor":INK,"axes.labelcolor":INK,"xtick.color":INK,"ytick.color":INK,
                     "axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":1.4,"savefig.dpi":300,"savefig.bbox":"tight"})
os.makedirs("fig",exist_ok=True)
J=lambda f: json.load(open(f))
# ---------- Fig 1: structures ----------
cpk={"Mo":("#5B7FA6",0.62),"S":("#E3B505",0.48)}
def draw(ax,xyz,title):
    L=open(xyz).read().split("\n")[2:]; at=[(l.split()[0],np.array(list(map(float,l.split()[1:4])))) for l in L if l.strip()]
    P=np.array([p for _,p in at]); c=P.mean(0); P=P-c
    # view: rotate slightly for depth
    th=np.deg2rad(-60); R=np.array([[1,0,0],[0,np.cos(th),-np.sin(th)],[0,np.sin(th),np.cos(th)]]); Q=P@R.T
    order=np.argsort(Q[:,2])
    for i in range(len(at)):
        for j in range(i+1,len(at)):
            if {at[i][0],at[j][0]}=={"Mo","S"} and np.linalg.norm(P[i]-P[j])<2.7:
                ax.plot([Q[i,0],Q[j,0]],[Q[i,1],Q[j,1]],color="#9AA3B2",lw=5,zorder=1)
    for i in order:
        col,r=cpk[at[i][0]]; ax.add_patch(plt.Circle(Q[i,:2],r,color=col,ec=INK,lw=1.2,zorder=2+Q[i,2]))
    ax.set_aspect("equal"); ax.axis("off"); ax.autoscale_view(); m=np.abs(Q[:,:2]).max()+0.9; ax.set_xlim(-m,m); ax.set_ylim(-m,m)
    ax.set_title(title,fontsize=20,color=INK,pad=6)
xyzs=[f for f in ["MoS2_opt.xyz","Mo3S6_opt.xyz"] if os.path.exists(f)]
fig,axs=plt.subplots(1,len(xyzs),figsize=(4.2*len(xyzs),4.2))
for ax,f in zip(np.atleast_1d(axs),xyzs): draw(ax,f,{"MoS2_opt.xyz":"MoS$_2$","Mo3S6_opt.xyz":"Mo$_3$S$_6$"}[f])
fig.savefig("fig/fig1_structures.png",transparent=True); plt.close(fig)
# ---------- Fig 2: correlation diagnostics ----------
sysl=[s for s in ["MoS2","Mo3S6"] if os.path.exists(f"{s}_classical.json")]
fig,(a1,a2)=plt.subplots(1,2,figsize=(12,4.6),gridspec_kw={"width_ratios":[1.15,1]})
w=0.38
for k,s in enumerate(sysl):
    d=J(f"{s}_classical.json"); occ=d["CAS44"]["nat_occ"]
    x=np.arange(len(occ))+(k-0.5*(len(sysl)-1))*w
    a1.bar(x,occ,w,color=[TEAL,GOLD][k],label={"MoS2":"MoS$_2$","Mo3S6":"Mo$_3$S$_6$"}[s])
a1.axhline(2,color=GREY,ls=":",lw=1.2); a1.axhline(0,color=GREY,lw=1)
a1.set_xticks(range(4)); a1.set_xticklabels([f"NO{i+1}" for i in range(4)]); a1.set_ylabel("Occupation, CAS(4e,4o)"); a1.legend(frameon=False)
a1.set_title("Natural-orbital occupations",fontsize=17,loc="left")
fx=["pbe","pbe0","hse06"]; 
for k,s in enumerate(sysl):
    d=J(f"{s}_classical.json"); a2.bar(np.arange(3)+(k-0.5*(len(sysl)-1))*w,[d[f"gap_{f}_eV"] for f in fx],w,color=[TEAL,GOLD][k])
a2.set_xticks(range(3)); a2.set_xticklabels(["PBE","PBE0","HSE06"]); a2.set_ylabel("HOMO–LUMO gap (eV)")
a2.set_title("Functional dependence of the gap",fontsize=17,loc="left")
fig.tight_layout(); fig.savefig("fig/fig2_correlation.png",transparent=True); plt.close(fig)
# ---------- Fig 3: accuracy vs circuit cost (noiseless) ----------
fig,ax=plt.subplots(figsize=(8,5.2))
for k,s in enumerate(sysl):
    f=f"{s}_cas44_adapt.json"
    if not os.path.exists(f): continue
    st=J(f)["steps"]; ax.semilogy([x["k"] for x in st],[max(abs(x["err_mHa"]),1e-3) for x in st],"o-",color=[TEAL,GOLD][k],lw=2.5,ms=9,
                                   label={"MoS2":"MoS$_2$ (4e,4o)","Mo3S6":"Mo$_3$S$_6$ (4e,4o)"}[s])
ax.axhspan(1e-3,1.6,color=TEAL,alpha=0.08); ax.axhline(1.6,color=TEAL,ls="--",lw=1.5); ax.text(0.0,1.25,"chemical accuracy (1.6 mHa)",ha="left",va="top",fontsize=13,color=TEAL)
ax.set_xlabel("ADAPT-VQE excitations in circuit"); ax.set_ylabel("|E$_{VQE}$ − E$_{CASCI}$| (mHa)"); ax.set_ylim(1e-2,300); ax.legend(frameon=False,loc="upper right")
ax.set_title("Accuracy vs. circuit size (noiseless)",fontsize=17,loc="left")
fig.tight_layout(); fig.savefig("fig/fig3_adapt.png",transparent=True); plt.close(fig)
# ---------- Fig 4: error-mitigation ladder ----------
rows=[]
for f in ["MoS2_cas44_A2_FakeBrisbane_noisy.json","MoS2_cas44_A2_FakeMarrakesh_noisy.json","Mo3S6_cas44_A2_FakeMarrakesh_noisy.json"]:
    if not os.path.exists(f): continue
    d=J(f); e=d["VQE"]["err_mHa_mean_std"]; bk=f.split("_A2_Fake")[1].split("_")[0]
    rows.append((d["system"],{"Brisbane":"Eagle","Marrakesh":"Heron"}[bk],e,d["VQE"]["err_ideal_mHa"]))
fig,ax=plt.subplots(figsize=(10,5.6))
labs=["raw","readout_mit","ZNE_linear"]; names=["Raw","+ TREX readout","+ ZNE (λ→0)"]; cols=[GREY,GOLD,TEAL]; W=0.8/3
for i,(s_,b,e,ideal) in enumerate(rows):
    for j,l in enumerate(labs):
        m,sd=e[l]; ax.bar(i+(j-1)*W,abs(m),W,yerr=sd,color=cols[j],capsize=5,label=names[j] if i==0 else None,error_kw={"lw":1.4})
        ax.text(i+(j-1)*W,(abs(m)+sd)*1.18,f"{m:.1f}",ha="center",fontsize=12,color=INK)
    ax.plot([i-0.45,i+0.45],[ideal,ideal],color=NAVY,lw=2.5,ls=(0,(4,2)),label="noiseless circuit" if i==0 else None)
ax.set_yscale("log"); ax.set_ylim(0.1,1000); ax.axhline(1.6,color=TEAL,ls=":",lw=1.5); ax.text(len(rows)-0.55,1.75,"1.6 mHa",color=TEAL,fontsize=12)
ax.set_xticks(range(len(rows))); ax.set_xticklabels([("MoS$_2$" if s_=="MoS2" else "Mo$_3$S$_6$")+f"\n{b} noise model" for s_,b,e,i_ in rows],fontsize=14)
ax.set_ylabel("|error| vs CASCI (mHa)"); ax.legend(frameon=False,ncol=4,loc="upper center",bbox_to_anchor=(0.5,1.13),fontsize=13)
fig.tight_layout(); fig.savefig("fig/fig4_mitigation.png",transparent=True); plt.close(fig)
# ---------- Fig 5: ZNE extrapolation (MoS2, Heron) ----------
f="MoS2_cas44_A2_FakeMarrakesh_noisy.json"
if os.path.exists(f):
    d=J(f); e=d["VQE"]["err_mHa_mean_std"]; lam=np.array([1,3,5]); y=np.array([e["readout_mit"][0],e["lam3_mit"][0],e["lam5_mit"][0]]); sd=[e["readout_mit"][1],e["lam3_mit"][1],e["lam5_mit"][1]]
    fig,ax=plt.subplots(figsize=(6.4,4.8)); p=np.polyfit(lam,y,1); xx=np.linspace(0,5.3,50)
    ax.plot(xx,np.polyval(p,xx),"--",color=TEAL,lw=2); ax.errorbar(lam,y,yerr=sd,fmt="o",color=NAVY,ms=10,capsize=5,label="measured (folded circuits)")
    ax.errorbar([0],[e["ZNE_linear"][0]],yerr=[e["ZNE_linear"][1]],fmt="*",color=TEAL,ms=20,capsize=5,label="zero-noise estimate")
    ax.axhspan(-1.6,1.6,color=TEAL,alpha=0.1); ax.axhline(0,color=GREY,lw=1)
    ax.set_xlabel("noise scale factor λ"); ax.set_ylabel("E − E$_{CASCI}$ (mHa)"); ax.legend(frameon=False,loc="upper left",fontsize=13)
    fig.tight_layout(); fig.savefig("fig/fig5_zne.png",transparent=True); plt.close(fig)
print(sorted(os.listdir("fig")))
