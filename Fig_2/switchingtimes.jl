using PyPlot, PyCall, MAT, LsqFit
gspec=pyimport("matplotlib.gridspec")
#########################################################################################
# Settings
#########################################################################################
# ----- Font style -----
rcParams=PyDict(matplotlib["rcParams"])
rcParams["text.usetex"]=true;
rcParams["axes.linewidth"]=1.25
PyPlot.rc("ytick", labelsize=12, left=true, right=true, direction="in")
PyPlot.rc("xtick", labelsize=12, top=true, bottom=true, direction="in")
PyPlot.rc("font", family="Helvetica", size=11)
PyPlot.rc("lines",linewidth=1.5,markersize=4)

#########################################################################################
# Critical slowing down - normal vs soliton (Figure 2b and 2c)
#########################################################################################
clr_m=(165,77,21)./255; clr_s=(41,120,120)./255
dir=joinpath(@__DIR__,"critical_loss_data/")

# Data files for the two gamma levels (G1, G3), each with a mirror (_m) and soliton (_s) run.
n1_m="A0Dm_ClosingGapG1.mat"; n1_s="A0D_ClosingGapG1.mat"
n3_m="A0Dm_ClosingGapG3.mat"; n3_s="A0D_ClosingGapG3.mat"
names=[n1_m,n1_s,n3_m,n3_s]
xpl=Vector{Vector{Float64}}(undef,4); ypl=Vector{Vector{Float64}}(undef,4); yErr=Vector{Matrix{Float64}}(undef,4)

for n ∈ eachindex(names)
    # load data
    vars=matread(dir*names[n])
    λu=vars["lambdaUP"]; λd=vars["lambdaDOWN"]; Jvec=vars["Jvec"]
    Errd=vars["ErrDOWN"]; Erru=vars["ErrUP"]; γ=vars["gamma"]
    ypl[n]=(λu.+λd)./γ; xpl[n]=Jvec./γ; yErr[n]=(Errd.+Erru)./γ
end

# ----- Article figure -----
rcParams["text.usetex"]=true;
rcParams["axes.linewidth"]=0.75
PyPlot.rc("ytick", labelsize=8, left=true, right=true, direction="in")
PyPlot.rc("xtick", labelsize=8, top=true, bottom=true, direction="in")
PyPlot.rc("font", family="Helvetica", size=11)
PyPlot.rc("lines",linewidth=1,markersize=3)
cm=1/2.54
spec=gspec.GridSpec(nrows=1,ncols=2,wspace=0.2)
figG=figure(figsize=(2.1*8.6/3*cm,1.2*8.6/3*cm))

# ----- Double-exponential fit shared by both panels -----
JvecG=Array(0.1:0.01:1.4)
@. modelL(x,p) = p[1]*exp(p[2]*x)+p[3]*exp(p[4]*x)

# ----- Panel (b): γ = 1 -----
axG1=figG.add_subplot(spec[1])
#axG1.errorbar(xpl[1],ypl[1],yErr[1],capsize=4.0,marker="D",linestyle="None",color=clr_m)
#axG1.errorbar(xpl[2],ypl[2],yErr[2],capsize=4.0,marker="o",linestyle="None",color=clr_s)
axG1.scatter(xpl[1],ypl[1],marker="D",color=clr_m)
axG1.scatter(xpl[2],ypl[2],marker="o",color=clr_s)
axG1.set_yscale("log")

fit1_m=LsqFit.curve_fit(modelL,xpl[1][1:end],ypl[1][1:end],[10,1,10,1.0])
fit1_s=LsqFit.curve_fit(modelL,xpl[2][1:end-2],ypl[2][1:end-2],[10,1,10,1.0])
axG1.plot(JvecG,modelL(JvecG,fit1_m.param),"--",color=clr_m)
axG1.plot(JvecG,modelL(JvecG,fit1_s.param),"-.",color=clr_s)
axG1.set_ylim([0.0101,1.5]); axG1.set_xlim([0.2,1.25])
axG1.set_yticks([0.01,0.1,1])
axG1.set_yticklabels([L"$0.01$",L"$0.1$",L"$1$"])

# ----- Panel (c): γ = 3 -----
axG2=figG.add_subplot(spec[2])
#axG2.errorbar(xpl[3],ypl[3],yErr[3],capsize=4.0,marker="D",linestyle="None",color=clr_m)
#axG2.errorbar(xpl[4],ypl[4],yErr[4],capsize=4.0,marker="o",linestyle="None",color=clr_s)
axG2.scatter(xpl[3],ypl[3],marker="D",color=clr_m)
axG2.scatter(xpl[4],ypl[4],marker="o",color=clr_s)
axG2.set_yscale("log")

fit2_m=LsqFit.curve_fit(modelL,xpl[3][1:end],ypl[3][1:end],[10,1,10,1.0])
fit2_s=LsqFit.curve_fit(modelL,xpl[4][1:end-2],ypl[4][1:end-2],[10,1,10,1.0])
axG2.plot(JvecG,modelL(JvecG,fit2_m.param),"--",color=clr_m)
axG2.plot(JvecG,modelL(JvecG,fit2_s.param),"-.",color=clr_s)
axG2.set_ylim([0.0101,1.5]); axG2.set_xlim([0.2,0.7])
axG2.set_xticks([0.2,0.4,0.6]); axG2.set_yticklabels([])

# ----- Shared labels -----
axG1.set_xlabel(L"$J/\gamma$",fontsize=9,labelpad=0.5); axG2.set_xlabel(L"$J/\gamma$",fontsize=8,labelpad=0.5);
axG1.set_ylabel(L"$\lambda/\gamma$",fontsize=9,labelpad=-3.0)
axG1.text(0.14, 0.09, L"$\rm{(b)}$", fontsize=10, horizontalalignment="center", verticalalignment="center", transform=axG1.transAxes)
axG2.text(0.14, 0.09, L"$\rm{(c)}$", fontsize=10, horizontalalignment="center", verticalalignment="center", transform=axG2.transAxes)
display(figG); figG.savefig("ClosingGap.pdf",format="pdf",bbox_inches="tight")
