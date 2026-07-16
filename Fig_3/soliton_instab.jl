using PyPlot, PyCall, DelimitedFiles
gspec=pyimport("matplotlib.gridspec")
pyimport("matplotlib.colors")
#########################################################################################
# This script produces three figures from Fig_3 data:
#   1) "Formation & decay of the phase profile"    -> Figure 3(a) and 3(b)
#   2) "Oscillatory instability during relaxation" -> Figure 3(c) and 3(d)
#   3) "Relaxation times vs damping"                -> Figure B2 (Appendix B)
#########################################################################################
# Settings
#########################################################################################
# ----- Font style -----
rcParams=PyDict(matplotlib["rcParams"])
rcParams["text.usetex"]=true;
rcParams["axes.linewidth"]=0.75
PyPlot.rc("ytick", labelsize=8, left=true, right=true, direction="in")
PyPlot.rc("xtick", labelsize=8, top=true, bottom=true, direction="in")
PyPlot.rc("font", family="Helvetica", size=10)
PyPlot.rc("lines",linewidth=1.0,markersize=3)
cm=1/2.54

#########################################################################################
# Formation & decay of the phase profile (Figure 3a and 3b)
#########################################################################################
dir1=joinpath(@__DIR__,"data_soliton_dynamics/") # shared data root, reused by all three figures below
custom_cmap = matplotlib.colors.LinearSegmentedColormap.from_list("custom",[(88, 0, 0)./255, (165, 77, 21)./255, (237, 197, 141)./255, (255, 255, 224)./255, (185, 214, 199)./255, (41, 120, 120)./255, (0, 50, 51)./255])

spec=gspec.GridSpec(nrows=2,ncols=1,hspace=0.11)
fig=figure(figsize=(0.46*8.6*cm,0.775*8.6*cm))

# ----- Panel (a): formation -----
dir2=joinpath(@__DIR__,"soliton_attractor/")
Data=Vector{Matrix{Float64}}(undef,3)
lbls=["phase","spacegrid_traj","timegrid"]
for i ∈ eachindex(lbls)
    f = readdlm(dir2*"Fig3c_"*lbls[i]*"10.txt")
    Data[i]=f
end

ax2=fig.add_subplot(spec[1])
pax=ax2.pcolormesh(Data[2][:,1].-51,Data[3][:,1],Data[1]',cmap=custom_cmap,vmin=-pi,vmax=pi,shading="gouraud")
ax2.set_ylabel(L"$\omega_r t$",labelpad=0.1)
ax2.text(0.05, 0.8, L"$\rm{(a)~formation}$", fontsize=10,  bbox=Dict(
        "facecolor" => "white",  # Background color
        "edgecolor" => "none",  # Border color (optional)
        "boxstyle"  => "round,pad=0.15",  # Rounded box with padding
        "alpha"     => 0.4              # Transparency
    ), transform=ax2.transAxes)
ax2.set_xlim([-40,40]); ax2.set_ylim([0,80]); ax2.set_xticklabels([])
ax2.set_yticks([0,30,60])

# ----- Panel (b): decay -----
Data=Vector{Matrix{Float64}}(undef,3)
lbls=["phase","spacegrid","timegrid"]
for i ∈ eachindex(lbls)
    f = readdlm(dir1*"Fig3c_"*lbls[i]*"1.txt")
    Data[i]=f
end

ax1=fig.add_subplot(spec[2])
pax=ax1.pcolormesh(Data[2][:,1].-51,Data[3][:,1],Data[1]',cmap=custom_cmap,vmin=-pi,vmax=pi,shading="gouraud")
ax1.plot([-40,40],[7.05,7.05],"k:",linewidth=1.25)
ax1.set_ylabel(L"$\omega_r t$",labelpad=0.15)
ax1.set_xlabel(L"$\mathrm{site}~j$",labelpad=0.15)
ax1.set_xlim([-40,40])
ax1.text(0.05, 0.8, L"$\rm{(b)~decay}$", fontsize=10,  bbox=Dict(
        "facecolor" => "white",  # Background color
        "edgecolor" => "none",  # Border color (optional)
        "boxstyle"  => "round,pad=0.15",  # Rounded box with padding
        "alpha"     => 0.4              # Transparency
    ), transform=ax1.transAxes)

# ----- Shared colorbar -----
cbar=colorbar(pax,ax=[ax1,ax2],ticks=[-pi,0,pi],drawedges=false,orientation="horizontal",location="top",pad=0.025)
cbar.ax.set_xticklabels([L"$-\pi$",L"$0$",L"$\pi$"],fontsize=8)
cbar.set_label(L"$\Phi_j - \overline{\Phi}$",labelpad=1,fontsize=10)
cbar.outline.set_visible(false)
display(fig); fig.savefig(joinpath(@__DIR__,"Spatial_Phase_Profile.pdf"),format="pdf",bbox_inches="tight")

#########################################################################################
# Oscillatory instability during relaxation (Figure 3c and 3d)
#########################################################################################
Data_a=Vector{Matrix{Float64}}(undef,4); clr1=(148, 60, 14)./255; clr2=(41,120,120)./255
Data_b=Vector{Matrix{Float64}}(undef,4)
lbls_a=["Phase_times","Phase","Density_times","Density"]
lbls_b=["j=1_Current_times","j=1_Current","j=-1_Current_times","j=-1_Current"]
for i ∈ eachindex(lbls_a)
    f = readdlm(dir1*"Fig3a_"*lbls_a[i]*".txt")
    Data_a[i]=f
    g = readdlm(dir1*"Fig3b_"*lbls_b[i]*".txt")
    Data_b[i]=g
end

PyPlot.rc("lines",linewidth=0.9,markersize=3)
spec=gspec.GridSpec(nrows=2,ncols=1,hspace=0.125)
fig=figure(figsize=(0.54*8.6*cm,0.6*8.6*cm))

# ----- Panel (c): phase & density -----
axR1=fig.add_subplot(spec[1]); axR1_b=axR1.twinx()
axR1.plot(Data_a[1],Data_a[2],color=clr1)
axR1_b.plot(Data_a[3],Data_a[4],color=clr2)
axR1_b.plot([7.05,7.05],[-0.2,1.4],"k:")
axR1.set_xlim([0,15]); axR1.set_xticklabels([])
axR1.set_ylim([0,pi]); axR1.set_yticks([0,pi]); axR1.set_yticklabels([L"$0$",L"$\pi$"])
axR1_b.set_ylim([0,1.25]); axR1_b.set_yticks([0,1]); axR1_b.set_yticklabels([L"$0$",L"$1$"])
axR1.tick_params(axis="y",labelcolor=clr1); axR1_b.tick_params(axis="y",labelcolor=clr2);
axR1.set_ylabel(L"$|\Delta\Phi_0(t)|$",fontsize=10,color=clr1); axR1_b.set_ylabel(L"N_0/N_f",fontsize=10,color=clr2);
axR1_b.yaxis.set_label_coords(1.05,0.44)
axR1.text(0.095, 0.8, L"$\rm{(c)}$", fontsize=10, horizontalalignment="center", verticalalignment="center", transform=axR1.transAxes)

# ----- Panel (d): site current -----
axR2=fig.add_subplot(spec[2])
axR2.plot(Data_b[1],Data_b[2],"k",label=L"$j\!=\!-1$")
axR2.plot(Data_b[3],Data_b[4],"k--",label=L"$j\!=\!1$")
axR2.plot([7.05,7.05],[-0.6,0.6],"k:")
axR2.set_ylim([-0.55,0.5]); axR2.set_yticks([-0.5,0,0.5]); axR2.set_yticklabels([L"$-0.5$",L"$0$",L"$0.5$"])
axR2.set_xlim([0,15])
axR2.set_xlabel(L"$\omega_r t$",fontsize=10,labelpad=-0.5); axR2.set_ylabel(L"$I_{j,0}/N_f$"); axR2.yaxis.set_label_coords(-0.1,0.5)
axR2.text(0.095, 0.12, L"$\rm{(d)}$", fontsize=10, horizontalalignment="center", verticalalignment="center", transform=axR2.transAxes)
axR2.legend(ncol=2,frameon=true,handlelength=1.0,fontsize=7.5,handletextpad=0.2,
    columnspacing=0.5,loc=(0.35,0.04),edgecolor="white",borderpad=0.2)

display(fig); fig.savefig(joinpath(@__DIR__,"Relaxation.pdf"),format="pdf",bbox_inches="tight")

#########################################################################################
# Relaxation times vs damping (Figure B2, Appendix B)
#########################################################################################
lbls=["xdata","ydata","yerr"]; Jv=[0.4,0.5,0.6]; mrkrs=["o","^","s"]
clrs=[(88,0,0)./255,(165,77,21)./255,(237,197,141)./255]
PyPlot.rc("ytick", labelsize=9, left=true, right=true, direction="in")
PyPlot.rc("xtick", labelsize=9, top=true, bottom=true, direction="in")
fig,ax=subplots(figsize=(0.725*8.6cm,0.55*8.6cm))
for j in eachindex(Jv)
    data=Vector{Vector{Float64}}(undef,3)
    for i in eachindex(lbls)
        f = open(dir1*"Fig3d_J="*string(Jv[j])*"_"*lbls[i]*".txt") do f
            readlines(f) |> (s->parse.(Float64, s))
        end
        data[i]=f
    end
    ax.errorbar(data[1],data[2],yerr=data[3],marker=mrkrs[j],fmt=" ",
        capsize=3,color=clrs[j])
end
ax.set_yscale("log"); ax.set_ylim([0.35,80]); ax.set_xlim([0,0.7])
ax.set_xlabel(L"$\gamma/J$"); ax.set_ylabel(L"$\Delta t_r(\omega_r^{-1})$")
ax.legend([L"$0.4$",L"$0.5$",L"$0.6$"],title=L"$J/\hbar\omega_r$",ncol=2,
    handletextpad=0.05,columnspacing=0.5,frameon=false,fontsize=9)
ax.plot([0,0.85],10 .^(2.7*[0,0.85]).- 0.45,"--",color=clrs[3])
display(fig); fig.savefig(joinpath(@__DIR__,"RelaxationTime.pdf"),format="pdf",bbox_inches="tight")
