using PyCall, PyPlot, MAT
gspec=pyimport("matplotlib.gridspec")
#########################################################################################
# Settings
#########################################################################################
# ----- Font style -----
rcParams=PyDict(matplotlib["rcParams"])
rcParams["text.usetex"]=true;
rcParams["axes.linewidth"]=1.0
PyPlot.rc("ytick", labelsize=12, left=true, right=true, direction="in")
PyPlot.rc("xtick", labelsize=12, top=true, bottom=true, direction="in")
PyPlot.rc("font", family="Helvetica", size=11)
PyPlot.rc("lines",linewidth=1.5,markersize=4)

#########################################################################################
# Stabilisation through losses (Figure B1, Appendix B)
#########################################################################################
dir=joinpath(@__DIR__,"critical_loss_data/")
varA=matread(dir*"A0D_CriticalLoss2.mat"); varT=matread(dir*"T0D_CriticalLoss3.mat")
ΩA=varA["Omega"]; GvecA=varA["Gvec"]; JvecA=varA["Jvec"]; U=varA["U"]
ΩT=varT["Omega"]; GvecT=varT["Gvec"]

rcParams["axes.linewidth"]=1.0
PyPlot.rc("ytick", labelsize=9, left=true, right=true, direction="in")
PyPlot.rc("xtick", labelsize=9, top=true, bottom=true, direction="in")
PyPlot.rc("font", family="Helvetica", size=10)
PyPlot.rc("lines",linewidth=1.25,markersize=4)

style=["-.","--",":","-"]
clr1=(165,77,21)./255; clr2=(41,120,120)./255; clr=[clr1,clr2]
cm=1/2.54
figγ,axγ=subplots(figsize=(0.7*8.6cm,0.55*8.6cm))
axγ.plot(GvecA/JvecA[2],imag.(ΩA[:,2,1])/U,style[1],color=clr[1],alpha=0.75,zorder=2)
axγ.plot(GvecA/JvecA[3],imag.(ΩA[:,3,1])/U,style[1],color=clr[2],alpha=0.75,zorder=3)
for n ∈ 2:3, j ∈ 2:3
    axγ.plot(GvecA/JvecA[j],imag.(ΩA[:,j,n])/U,style[n],color=clr[j-1],alpha=0.75,zorder=j)
end
axγ.plot(GvecT/JvecA[2],imag.(ΩT[:,2,3])/U,style[end],color=clr[1],zorder=4,label=L"$J/\mu=0.14$")
axγ.plot(GvecT/JvecA[3],imag.(ΩT[:,3,3])/U,style[end],color=clr[2],zorder=4,label=L"$J/\mu=0.34$")

axγ.set_yscale("log"); axγ.set_ylim([1e-4,5]);
axγ.set_yticks([10^(-4),0.01,1])
axγ.set_ylabel(L"$\mathrm{max}(\mathrm{Im}[\omega_l])/U$"); axγ.set_xlabel(L"$\gamma/J$",labelpad=0.2)
axγ.set_xlim([0.2,13]); axγ.set_xscale("log")

axγ.legend(frameon=false,loc="upper right",handlelength=1.1,handletextpad=0.4)
display(figγ); figγ.savefig(joinpath(@__DIR__,"CriticalLoss.pdf"),format="pdf",bbox_inches="tight")