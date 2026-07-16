include("source.jl")
using PyPlot, PyCall, MAT
gspec=pyimport("matplotlib.gridspec")
#########################################################################################
# Settings
#########################################################################################
# ----- Figure size -----
standard=(3.38,3.5) # width and heigt in inches
wdth=1.85; hght=1.90

# ----- Fontsize -----
groot=16; mid=14; klein=12;

# ----- Font style -----
rcParams=PyDict(matplotlib["rcParams"])
rcParams["text.usetex"]=true;
rcParams["axes.linewidth"]=1.25
PyPlot.rc("ytick", labelsize=12, left=true, right=true, direction="in")
PyPlot.rc("xtick", labelsize=12, top=true, bottom=true, direction="in")
PyPlot.rc("font", family="Helvetica", size=11)
PyPlot.rc("lines",linewidth=1.5,markersize=4)

# ----- Colors -----
#bluegreen
grgr=(112/255,177/255,130/255,1)
lgr=(171/255,226/255,199/255,1)
grbl=(237/255,253/255,214/255,1)
lbl=(145/255,201/255,214/255,1)
blbl=(87/255,156/255,173/255,1)
#...

darkblue=(31,102,169)./255
medblue=(52,148,204)./255
lightblue=(141,197,228)./255

darkteal=(31,111,112)./255
medteal=(84,161,162)./255
lightteal=(159,200,200)./255

darkred=(159,0,0)./255
medred=(196,103,102)./255
lightred=(216,165,166)./255

darkyel=(219,161,28)./255

black=(0,0,0)
lightgray=(0.6,0.6,0.6)
darkgray=(0.25,0.25,0.25)


#########################################################################################
# Compare PhaseDiagram
#########################################################################################
dir="D:\\Storage\\Driv_Disp_BEC\\Data\\3D\\WN_exphyst\\"
dirM=dir*"PhaseDiagram_mirror\\"
dirS=dir*"PhaseDiagram_sol\\"
dirs=[dirS;dirM]; Jdim=40
nms=["A2D_RelaxG";"A2Dm_RelaxG"]; Glevel=[1,3,5,7,10,15,20,25,30,40,50]; ΔN=zeros(Float64,Jdim,size(Glevel,1),2)
γvec=zeros(Float64,size(Glevel,1),2); Jvec=zeros(Float64,Jdim,2); Nf=ones(Float64,Jdim,2)
for cnf = 1:2, g ∈ eachindex(Glevel)
    nameG=nms[cnf]*string(Glevel[g])*"J";
    for j = 2:Jdim
        varsE=matread(dirs[cnf]*nameG*string(j)*"_E.mat");
        varsF=matread(dirs[cnf]*nameG*string(j)*"_F.mat");
        Nₓ=varsE["Nx"]; Nʸ=varsE["Ny"]; Sx=Nₓ+1; Sy=Nʸ+1; Pcut=varsE["Pcut"]
        nn=varsE["nn"]; cnf==2 ? (cz=0) : (cz=1)
 
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
#cb.set_label(L"$\Delta N/N_f$",labelpad=5)
#cb.ax.set_ylabel(L"$\Delta N/N_f$",fontsize=8,rotation=0,labelpad=8); #cb.ax.yaxis.set_label_coords(0.1,0.1)
# shift colorbar slightly to the left
pos = cb.ax.get_position()
cb.ax.set_position([pos.x0-0.075, pos.y0, pos.width, pos.height])

# put label to the right of the bar
cb.ax.set_xlabel(L"$\Delta N/N_f$", fontsize=10)
cb.ax.xaxis.set_label_coords(1.13,0.25)

varsEx=matread("ExpData.mat")
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
#ax2.text(0.05,0.4, L"\textbf{superfluid}",rotation=25,color="w",fontsize=14)
ax1.set_xticks([0,0.2,0.4]); ax2.set_xticks([0,0.2,0.4])

# ----- Transition points -----
#sing CurveFit
#indT=zeros(Int64,size(Glevel))
#for g ∈ eachindex(Glevel)
#    indT[g]=findfirst(x->x>0.2,ΔN[:,g,2]./Nf[:,2])
#end
#a,b=linear_fit(γvec[:,2],Jvec[indT,2])
#ax2.plot(γvec[:,2],2/3*γvec[:,2],"-o",color="w")

#ax1.set_xlabel(L"$\gamma/\omega_r$");
ax1.set_xlim([0,0.55]); ax1.set_ylim([0,0.55])
ax2.set_xlim([0,0.55]); ax2.set_ylim([0,0.55])
ax2.set_xlabel(L"$\gamma/\hbar\omega_r$",fontsize=10,labelpad=0.5); #ax1.set_xticklabels([])
ax1.set_ylabel(L"$J/\hbar\omega_r$",fontsize=10); ax2.set_yticklabels([])
ax1.set_xlabel(L"$\gamma/\hbar\omega_r$",fontsize=10,labelpad=0.5);
ax1.text(0.12, 0.92, L"$\rm{(a)}$", fontsize=10, horizontalalignment="center", verticalalignment="center", color="white",transform=ax1.transAxes)
ax2.text(0.12, 0.92, L"$\rm{(b)}$", fontsize=10, horizontalalignment="center", verticalalignment="center", color="white",transform=ax2.transAxes)
display(fig); fig.savefig("PhaseDiagram3D_b.pdf",format="pdf",bbox_inches="tight")


#########################################################################################
# Compare PhaseDiagram: Supplemental
#########################################################################################
dir="D:\\Storage\\Driv_Disp_BEC\\Data\\3D\\WN_exphyst\\"
dir1=dir*"PhaseDiagram_sol\\"; dir2=dir*"PhaseDiagram_L70\\";
nms=[dir2*"A2D_N1p5_RelaxG";dir2*"A2D_N2_RelaxG"];
Jdim=30; Glevel=[1,3,5,7,10,12,15,17,20]; ΔN=zeros(Float64,Jdim,size(Glevel,1),3)
γvec=zeros(Float64,size(Glevel,1),3); Jvec=zeros(Float64,Jdim,3); Nf=ones(Float64,Jdim,3)
for cnf = 1:2, g ∈ eachindex(Glevel)
    nameG=nms[cnf]*string(Glevel[g])*"J";
    for j = 1:Jdim
        varsE=matread(nameG*string(j)*"_E.mat");
        varsF=matread(nameG*string(j)*"_F.mat");
        Nₓ=varsE["Nx"]; Nʸ=varsE["Ny"]; Sx=Nₓ+1; Sy=Nʸ+1; Pcut=varsE["Pcut"]
        nn=varsE["nn"]; cz=1
 
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
#cb.set_label(L"$\Delta N/N_f$",labelpad=5)
#cb.ax.set_ylabel(L"$\Delta N/N_f$",fontsize=8,rotation=0,labelpad=8); #cb.ax.yaxis.set_label_coords(0.1,0.1)
# shift colorbar slightly to the left
pos = cb.ax.get_position()
cb.ax.set_position([pos.x0-0.075, pos.y0, pos.width, pos.height])

# put label to the right of the bar
cb.ax.set_xlabel(L"$\Delta N/N_f$", fontsize=10)
cb.ax.xaxis.set_label_coords(1.13,0.25)

varsEx=matread("ExpData.mat")
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

# ----- Transition points -----
#sing CurveFit
#indT=zeros(Int64,size(Glevel))
#for g ∈ eachindex(Glevel)
#    indT[g]=findfirst(x->x>0.2,ΔN[:,g,2]./Nf[:,2])
#end
#a,b=linear_fit(γvec[:,2],Jvec[indT,2])
#ax2.plot(γvec[:,2],2/3*γvec[:,2],"-o",color="w")

#ax1.set_xlabel(L"$\gamma/\omega_r$");
ax1.set_xlim([0,0.215]); ax1.set_ylim([0,0.375])
ax2.set_xlim([0,0.215]); ax2.set_ylim([0,0.375])
ax2.set_xlabel(L"$\gamma/\hbar\omega_r$",fontsize=10,labelpad=0.5); #ax1.set_xticklabels([])
ax1.set_ylabel(L"$J/\hbar\omega_r$",fontsize=10); ax2.set_yticklabels([])
ax1.set_xlabel(L"$\gamma/\hbar\omega_r$",fontsize=10,labelpad=0.5);
ax1.text(0.12, 0.92, L"$\rm{(a)}$", fontsize=10, horizontalalignment="center", verticalalignment="center", color="white",transform=ax1.transAxes)
ax2.text(0.12, 0.92, L"$\rm{(b)}$", fontsize=10, horizontalalignment="center", verticalalignment="center", color="white",transform=ax2.transAxes)
display(fig); fig.savefig("PhaseDiagram_L70.pdf",format="pdf",bbox_inches="tight")


#########################################################################################
# PhaseDiagram
#########################################################################################
#dir="C:\\Users\\robbe\\Documents\\Driven_Disp_BEC\\3Dimensional\\output_files\\"
dir="D:\\Storage\\Driv_Disp_BEC\\Data\\3D\\WN_exphyst\\PhaseDiagram_L70\\"
nms="A2D_N2_RelaxG"; n_J=30; Glevel=[1,3,5,7,10,12,15,17,20]; ΔN=zeros(Float64,n_J,size(Glevel,1))
γvec=zeros(Float64,size(Glevel,1)); Jvec=zeros(Float64,n_J); Nf=ones(Float64,n_J)
for g ∈ eachindex(Glevel)
    nameG=nms*string(Glevel[g])*"J";
    for j = 1:n_J
        varsE=matread(dir*nameG*string(j)*"_E.mat");
        varsF=matread(dir*nameG*string(j)*"_F.mat");
        Nₓ=varsE["Nx"]; Nʸ=varsE["Ny"]; Sx=Nₓ+1; Sy=Nʸ+1; Pcut=varsE["Pcut"]
        nn=varsE["nn"]; cz=1
 
        Cend_E=varsE["Ct"][cz*Sx+1:(cz+1)*Sx,:,end,:]; nsim_E=size(Cend_E,3);
        Cend_F=varsF["Ct"][cz*Sx+1:(cz+1)*Sx,:,end,:]; nsim_F=size(Cend_F,3); 
        j==n_J && (γvec[g]=varsE["gamma"])
        g==1 && (Jvec[j]=varsE["J"])
        Cn2_E=dropdims(sum(abs2.(Cend_E),dims=3),dims=3)/nsim_E.-0.5*Pcut[1:Sx,:]
        Cn2_F=dropdims(sum(abs2.(Cend_F),dims=3),dims=3)/nsim_F.-0.5*Pcut[1:Sx,:]
        ΔN[j,g]=abs(sum(Cn2_F)-sum(Cn2_E));

        if g==1
            Cstrt_F=varsF["Ct"][cz*Sx+1:(cz+1)*Sx,:,1,:]
            Cn2_F=dropdims(sum(abs2.(Cstrt_F),dims=3),dims=3)/nsim_F.-0.5*Pcut[1:Sx,:]
            Nf[j]=sum(Cn2_F)
        end
    end
end
# ----- Article figure -----
pyimport("matplotlib.colors")
custom_cmap = matplotlib.colors.LinearSegmentedColormap.from_list("custom",[(88, 0, 0)./255, (165, 77, 21)./255, (237, 197, 141)./255,
    (255, 255, 224)./255, (185, 214, 199)./255, (41, 120, 120)./255, (0, 50, 51)./255])

fig=figure(figsize=(2.5wdth,3hght)); ax1=fig.add_subplot(); 
pc1=ax1.pcolormesh(γvec[:],Jvec[:],ΔN[:,:]./Nf[:],vmin=0,vmax=0.8,edgecolors="face",shading="nearest")
ax1.set_xlim([0,0.55]); ax1.set_ylim([0,0.55])
cb=fig.colorbar(pc1,ax=ax1,
        orientation="vertical",location="right",aspect=30,shrink=0.75)
cb.set_label(label=L"$\Delta N/N_f$",size=15)
cb.outline.set_visible(false)
ax1.plot([0,0.5],[0,0.5/4],"--",color="white")
#display(fig)

varsEx=matread("ExpData.mat")
Xlb=varsEx["Xlb"][:]; Ylb=varsEx["Ylb"][:]; ΔXlb=varsEx["DeltaXlb"][:]; ΔYlb=varsEx["DeltaYlb"][:]
Xub=varsEx["Xub"][:]; Yub=varsEx["Yub"][:]; ΔXub=varsEx["DeltaXub"][:]; ΔYub=varsEx["DeltaYub"][:]
ax1.errorbar(Xub,Yub,ΔYub,ΔXub,marker="s",linestyle="none",ecolor=darkred,mec=darkred,mfc=lightred)
ax1.errorbar(Xlb,Ylb,ΔYlb,ΔXlb,marker="o",linestyle="none",ecolor=darkblue,mec=darkblue,mfc=lightblue)
Xftub=varsEx["Xftub"][:]; Yftub=varsEx["Yftub"][:]
#ax1.plot(Xftub,Yftub,color=darkred);
xJ=[0.0,0.6]; γJ=1.0*xJ; γJ4=xJ/4;
ax1.plot(xJ,γJ,"--",color=darkblue);
#ax1.plot(xJ,γJ4,"-.",color="w"); ax2.plot(xJ,γJ4,"-.",color="w")
#ax2.text(0.3, 0.032, L"$J=\hbar\gamma/4$", color="w",rotation=15,rotation_mode="anchor")
#ax2.text(0.05,0.4, L"\textbf{superfluid}",rotation=25,color="w",fontsize=14)

# ----- Transition points -----
using CurveFit
indT=zeros(Int64,size(Glevel))
for g ∈ eachindex(Glevel)
    indT[g]=findfirst(x->x>0.1,ΔN[:,g]./Nf[:])
end
a,b=linear_fit(γvec[:],Jvec[indT])
#ax1.plot(Array(0:0.01:0.5),a.+b*Array(0:0.01:0.5),"-",color="w")

ax1.set_xlim([0,0.25]); ax1.set_ylim([0,0.475])
ax1.set_xlabel(L"$\gamma/\omega_r$"); ax1.set_ylabel(L"$J/\hbar\omega_r$",fontsize=15);
#ax1.text(0.08, 0.93, "(a)", horizontalalignment="center", verticalalignment="center", color="white",transform=ax1.transAxes)
#ax2.text(0.08, 0.93, "(b)", horizontalalignment="center", verticalalignment="center", color="white",transform=ax2.transAxes)
ax1.set_title(L"$\sigma^2 = 4$")
display(fig);