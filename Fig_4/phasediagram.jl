using PyPlot, PyCall, MAT
gspec=pyimport("matplotlib.gridspec")
#########################################################################################
# This script produces two phase-diagram figures from TWA simulation data:
#   1) "Compare PhaseDiagram"              -> Figure 4 (main text)
#   2) "Compare PhaseDiagram: Supplemental" -> Figure E1 (Appendix E)
# Both blocks follow the same pattern: load the TWA data for a set of (gamma, J)
# points, compute the population-difference order parameter, then plot it as a
# two-panel colormesh with the experimental phase boundary overlaid.
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

# ----- Colors -----
darkblue=(31,102,169)./255
medblue=(52,148,204)./255
lightblue=(141,197,228)./255

darkred=(159,0,0)./255
lightred=(216,165,166)./255


#########################################################################################
# Compare PhaseDiagrams
#########################################################################################
dirM=joinpath(@__DIR__,"twa_data/PhaseDiagram_mirror/")
dirS=joinpath(@__DIR__,"twa_data/PhaseDiagram_sol/")
dirs=[dirS;dirM]; Jdim=40
nms=["A2D_RelaxG";"A2Dm_RelaxG"]; Glevel=[1,3,5,7,10,15,20,25,30,40,50]; ΔN=zeros(Float64,Jdim,size(Glevel,1),2)
γvec=zeros(Float64,size(Glevel,1),2); Jvec=zeros(Float64,Jdim,2); Nf=ones(Float64,Jdim,2)
# For each config (soliton/mirror), G-level and J, load the _E/_F datasets and
# compute the population difference ΔN between them, normalized by the initial
# population Nf.
for cnf = 1:2, g ∈ eachindex(Glevel)
    nameG=nms[cnf]*string(Glevel[g])*"J";
    for j = 2:Jdim
        varsE=matread(dirs[cnf]*nameG*string(j)*"_E.mat");
        varsF=matread(dirs[cnf]*nameG*string(j)*"_F.mat");
        Nₓ=varsE["Nx"]; Sx=Nₓ+1; Pcut=varsE["Pcut"]
        cnf==2 ? (cz=0) : (cz=1)

        Cend_E=varsE["Ct"][cz*Sx+1:(cz+1)*Sx,:,end,:]; nsim_E=size(Cend_E,3);
        Cend_F=varsF["Ct"][cz*Sx+1:(cz+1)*Sx,:,end,:]; nsim_F=size(Cend_F,3); 
        j==30 && (γvec[g,cnf]=varsE["gamma"])
        g==1 && (Jvec[j,cnf]=varsE["J"])
        Cn2_E=dropdims(sum(abs2.(Cend_E),dims=3),dims=3)/nsim_E.-0.5*Pcut[1:Sx,:]
        Cn2_F=dropdims(sum(abs2.(Cend_F),dims=3),dims=3)/nsim_F.-0.5*Pcut[1:Sx,:]
        ΔN[j,g,cnf]=abs(sum(Cn2_F)-sum(Cn2_E));

        if g==1
            Cstrt_F=varsF["Ct"][cz*Sx+1:(cz+1)*Sx,:,1,:]
            Cn2_F=dropdims(sum(abs2.(Cstrt_F),dims=3),dims=3)/nsim_F.-0.5*Pcut[1:Sx,:]
            Nf[j,cnf]=sum(Cn2_F)
        end
    end
end
# ----- Article figure -----
rcParams["axes.linewidth"]=0.75
PyPlot.rc("ytick", labelsize=8, left=true, right=true, direction="in")
PyPlot.rc("xtick", labelsize=8, top=true, bottom=true, direction="in")
PyPlot.rc("font", family="Helvetica", size=10)
PyPlot.rc("lines",linewidth=1.0,markersize=3.6)

cm=1/2.54
spec=gspec.GridSpec(ncols=2, nrows=1, wspace=0.08)
fig=figure(figsize=(1.1*8.6*cm,0.65*8.6*cm)); ax1=fig.add_subplot(spec[1]); ax2=fig.add_subplot(spec[2])
pc1=ax1.pcolormesh(γvec[:,1],Jvec[:,1],ΔN[:,:,1]./Nf[:,1],vmin=0,vmax=0.875,edgecolors="face",shading="nearest")
pc2=ax2.pcolormesh(γvec[:,2],Jvec[:,2],ΔN[:,:,2]./Nf[:,2],vmin=0,vmax=0.875,edgecolors="face",shading="nearest")

cb=fig.colorbar(pc1,ax=[ax1,ax2],
        orientation="horizontal",location="top",aspect=27.5,shrink=0.75)
cb.outline.set_visible(false)
# shift colorbar slightly to the left
pos = cb.ax.get_position()
cb.ax.set_position([pos.x0-0.075, pos.y0, pos.width, pos.height])

# put label to the right of the bar
cb.ax.set_xlabel(L"$\Delta N/N_f$", fontsize=10)
cb.ax.xaxis.set_label_coords(1.13,0.25)

varsEx=matread(joinpath(@__DIR__,"twa_data/ExpData.mat"))
Xlb=varsEx["Xlb"][:]; Ylb=varsEx["Ylb"][:]; ΔXlb=varsEx["DeltaXlb"][:]; ΔYlb=varsEx["DeltaYlb"][:]
Xub=varsEx["Xub"][:]; Yub=varsEx["Yub"][:]; ΔXub=varsEx["DeltaXub"][:]; ΔYub=varsEx["DeltaYub"][:]
ax1.errorbar(Xub,Yub,ΔYub,ΔXub,marker="s",linestyle="none",ecolor=darkred,mec=darkred,mfc=lightred)
ax1.errorbar(Xlb,Ylb,ΔYlb,ΔXlb,marker="o",linestyle="none",ecolor=darkblue,mec=darkblue,mfc=lightblue)
ax2.errorbar(Xub,Yub,ΔYub,ΔXub,marker="s",linestyle="none",ecolor=darkred,mec=darkred,mfc=lightred)
ax2.errorbar(Xlb,Ylb,ΔYlb,ΔXlb,marker="o",linestyle="none",ecolor=darkblue,mec=darkblue,mfc=lightblue)
Xftub=varsEx["Xftub"][:]; Yftub=varsEx["Yftub"][:]
ax1.plot(Xftub,Yftub,color=darkred); ax2.plot(Xftub,Yftub,color=darkred)
xJ=[0.0,0.6]; γJ=1.0*xJ; γJ4=xJ/4;
ax1.plot(xJ,γJ,"--",color=medblue); ax2.plot(xJ,γJ,"--",color=medblue);
ax1.plot(xJ,γJ4,"-.",color="w"); ax2.plot(xJ,γJ4,"-.",color="w")
ax2.text(0.3, 0.030, L"$J=\gamma/4$", fontsize=8, color="w",rotation=15,rotation_mode="anchor")
ax1.set_xticks([0,0.2,0.4]); ax2.set_xticks([0,0.2,0.4])

ax1.set_xlim([0,0.55]); ax1.set_ylim([0,0.55])
ax2.set_xlim([0,0.55]); ax2.set_ylim([0,0.55])
ax2.set_xlabel(L"$\gamma/\hbar\omega_r$",fontsize=10,labelpad=0.5)
ax1.set_ylabel(L"$J/\hbar\omega_r$",fontsize=10); ax2.set_yticklabels([])
ax1.set_xlabel(L"$\gamma/\hbar\omega_r$",fontsize=10,labelpad=0.5);
ax1.text(0.12, 0.92, L"$\rm{(a)}$", fontsize=10, horizontalalignment="center", verticalalignment="center", color="white",transform=ax1.transAxes)
ax2.text(0.12, 0.92, L"$\rm{(b)}$", fontsize=10, horizontalalignment="center", verticalalignment="center", color="white",transform=ax2.transAxes)
display(fig); fig.savefig("PhaseDiagram3D_b.pdf",format="pdf",bbox_inches="tight")


#########################################################################################
# Compare PhaseDiagrams: Supplemental (Figure E1, Appendix E)
#########################################################################################
dir2=joinpath(@__DIR__,"twa_data/PhaseDiagram_L70/")
nms=[dir2*"A2D_N1p5_RelaxG";dir2*"A2D_N2_RelaxG"];
Jdim=30; Glevel=[1,3,5,7,10,12,15,17,20]; ΔN=zeros(Float64,Jdim,size(Glevel,1),2)
γvec=zeros(Float64,size(Glevel,1),2); Jvec=zeros(Float64,Jdim,2); Nf=ones(Float64,Jdim,2)
# Same population-difference computation as above, for the two N-configurations (N=1.5, N=2).
for cnf = 1:2, g ∈ eachindex(Glevel)
    nameG=nms[cnf]*string(Glevel[g])*"J";
    for j = 1:Jdim
        varsE=matread(nameG*string(j)*"_E.mat");
        varsF=matread(nameG*string(j)*"_F.mat");
        Nₓ=varsE["Nx"]; Sx=Nₓ+1; Pcut=varsE["Pcut"]
        cz=1

        Cend_E=varsE["Ct"][cz*Sx+1:(cz+1)*Sx,:,end,:]; nsim_E=size(Cend_E,3);
        Cend_F=varsF["Ct"][cz*Sx+1:(cz+1)*Sx,:,end,:]; nsim_F=size(Cend_F,3); 
        j==30 && (γvec[g,cnf]=varsE["gamma"])
        g==1 && (Jvec[j,cnf]=varsE["J"])
        Cn2_E=dropdims(sum(abs2.(Cend_E),dims=3),dims=3)/nsim_E.-0.5*Pcut[1:Sx,:]
        Cn2_F=dropdims(sum(abs2.(Cend_F),dims=3),dims=3)/nsim_F.-0.5*Pcut[1:Sx,:]
        ΔN[j,g,cnf]=abs(sum(Cn2_F)-sum(Cn2_E));

        if g==1
            Cstrt_F=varsF["Ct"][cz*Sx+1:(cz+1)*Sx,:,1,:]
            Cn2_F=dropdims(sum(abs2.(Cstrt_F),dims=3),dims=3)/nsim_F.-0.5*Pcut[1:Sx,:]
            Nf[j,cnf]=sum(Cn2_F)
        end
    end
end
# ----- Article figure -----
rcParams["axes.linewidth"]=0.75
PyPlot.rc("ytick", labelsize=8, left=true, right=true, direction="in")
PyPlot.rc("xtick", labelsize=8, top=true, bottom=true, direction="in")
PyPlot.rc("font", family="Helvetica", size=10)
PyPlot.rc("lines",linewidth=1.0,markersize=3.6)

cm=1/2.54
spec=gspec.GridSpec(ncols=2, nrows=1, wspace=0.08)
fig=figure(figsize=(1.1*8.6*cm,0.65*8.6*cm)); ax1=fig.add_subplot(spec[1]); ax2=fig.add_subplot(spec[2])
pc1=ax1.pcolormesh(γvec[:,1],Jvec[:,1],ΔN[:,:,1]./Nf[:,1],vmin=0,vmax=0.875,edgecolors="face",shading="nearest")
pc2=ax2.pcolormesh(γvec[:,2],Jvec[:,2],ΔN[:,:,2]./Nf[:,2],vmin=0,vmax=0.875,edgecolors="face",shading="nearest")

cb=fig.colorbar(pc1,ax=[ax1,ax2],
        orientation="horizontal",location="top",aspect=27.5,shrink=0.75)
cb.outline.set_visible(false)
# shift colorbar slightly to the left
pos = cb.ax.get_position()
cb.ax.set_position([pos.x0-0.075, pos.y0, pos.width, pos.height])

# put label to the right of the bar
cb.ax.set_xlabel(L"$\Delta N/N_f$", fontsize=10)
cb.ax.xaxis.set_label_coords(1.13,0.25)

varsEx=matread(joinpath(@__DIR__,"twa_data/ExpData.mat"))
Xlb=varsEx["Xlb"][:]; Ylb=varsEx["Ylb"][:]; ΔXlb=varsEx["DeltaXlb"][:]; ΔYlb=varsEx["DeltaYlb"][:]
Xub=varsEx["Xub"][:]; Yub=varsEx["Yub"][:]; ΔXub=varsEx["DeltaXub"][:]; ΔYub=varsEx["DeltaYub"][:]
ax1.errorbar(Xub,Yub,ΔYub,ΔXub,marker="s",linestyle="none",ecolor=darkred,mec=darkred,mfc=lightred)
ax1.errorbar(Xlb,Ylb,ΔYlb,ΔXlb,marker="o",linestyle="none",ecolor=darkblue,mec=darkblue,mfc=lightblue)
ax2.errorbar(Xub,Yub,ΔYub,ΔXub,marker="s",linestyle="none",ecolor=darkred,mec=darkred,mfc=lightred)
ax2.errorbar(Xlb,Ylb,ΔYlb,ΔXlb,marker="o",linestyle="none",ecolor=darkblue,mec=darkblue,mfc=lightblue)
Xftub=varsEx["Xftub"][:]; Yftub=varsEx["Yftub"][:]
ax1.plot(Xftub,Yftub,color=darkred); ax2.plot(Xftub,Yftub,color=darkred)
xJ=[0.0,0.6]; γJ=1.0*xJ; γJ4=xJ/4;
ax1.plot(xJ,γJ,"--",color=medblue); ax2.plot(xJ,γJ,"--",color=medblue);
ax1.plot(xJ,γJ4,"-.",color="w"); ax2.plot(xJ,γJ4,"-.",color="w")
#ax2.text(0.3, 0.030, L"$J=\gamma/4$", fontsize=8, color="w",rotation=15,rotation_mode="anchor")
#ax2.text(0.05,0.4, L"\textbf{superfluid}",rotation=25,color="w",fontsize=14)
#ax1.set_xticks([0,0.2,0.4]); ax2.set_xticks([0,0.2,0.4])

ax1.set_xlim([0,0.215]); ax1.set_ylim([0,0.375])
ax2.set_xlim([0,0.215]); ax2.set_ylim([0,0.375])
ax2.set_xlabel(L"$\gamma/\hbar\omega_r$",fontsize=10,labelpad=0.5)
ax1.set_ylabel(L"$J/\hbar\omega_r$",fontsize=10); ax2.set_yticklabels([])
ax1.set_xlabel(L"$\gamma/\hbar\omega_r$",fontsize=10,labelpad=0.5);
ax1.text(0.12, 0.92, L"$\rm{(a)}$", fontsize=10, horizontalalignment="center", verticalalignment="center", color="white",transform=ax1.transAxes)
ax2.text(0.12, 0.92, L"$\rm{(b)}$", fontsize=10, horizontalalignment="center", verticalalignment="center", color="white",transform=ax2.transAxes)
display(fig); fig.savefig("PhaseDiagram_L70.pdf",format="pdf",bbox_inches="tight")