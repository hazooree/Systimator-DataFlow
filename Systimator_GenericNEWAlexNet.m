clear, close all
Line_Width = 1;Marker_Size = 30;Font_Size = 10;%pattern = ['o','s','d','v','p','h'];
CNN = [
    %Filters   Size    FM    IFM_Depth
    96       11       227      3
    256       5       27      96
    384       3       13     256
    384       3       13     384
    256       3      13     384];
HWR = [220 4.9];N_DSP = HWR(1);M_BRAM = 0.9*HWR(2);
F=16;O=6;P=4;Q=4;L=size(CNN,1);s=2;W=1;
Mb = 1/2^20;N_b = 16;
n_f=CNN(:,1);r_f=CNN(:,2);c_f=r_f;r=CNN(:,3);c=r;ch=CNN(:,4);
for rho = 0:1
    i = 0;
    for o=1:O
        for p=1:P
            for q=1:Q
                i = i + 1;
                c_sa(p) = 2.^p;ch_sa(q) = 2.^q;r_sa(q) = ch_sa(q)*max(r_f);
                n_dsp(i) = c_sa(p)*r_sa(q);
                A(i) = c_sa(p);
                B(i) = r_sa(q);
                for l=1:L
                    r_t(l,o) = min(ceil(r(1)./o/F),r(l));
                    c_t(l,o) = c(l);
                    d_H(l,o) = r_t(l,o)-r_f(l)+1;
                    d_V(l,o) = c_t(l,o)-c_f(l)+1;
                    M_FM(l,i)   = r_t(l,o)*c_t(l,o)*ch_sa(q)*N_b*Mb;
                    M_SP(l,i)   = c_sa(p)*N_b*Mb;
                    M_SA(l,i)   = ((1-rho)*n_f(l)+rho*c_sa(p))*d_V(l,o)*d_H(l,o)*N_b*Mb;
                    M_out(l,i) = M_SA(l,i)/s^2;
                    M_W_SA(l,i) = r_f(l)*(n_f(l)*(1-rho)+rho*c_sa(p))*ch_sa(q)*N_b*Mb;
                    Alpha(l,i,rho+1)  = ceil(n_f(l)/c_sa(p));
                    Beta(l,i)   = ceil(r(l)/r_t(l,o));
                    Gemma(l,i)  = ceil(ch(l)/ch_sa(q));
                    Omega(l,i)  = Alpha(l,i,rho+1)*Beta(l,i)*Gemma(l,i);
                    T_FM(l,i,rho+1)   = 1/W*(Alpha(l,i,rho+1)*rho+1-rho)*Beta(l,i)*Gemma(l,i)*M_FM(l,i)/Mb;
                    T_W(l,i,rho+1)    = 1/W*(Alpha(l,i,rho+1)*(1-rho)+rho)*Beta(l,i)*Gemma(l,i)*M_W_SA(l,i)/Mb;
                    T_SP(l,i)   = Omega(l,i)*(d_V(l,o)*d_H(l,o)+r_sa(q)-1)*r_f(l);
                    T_SA(l,i)   = T_SP(l, i)+Omega(l,i)*c_sa(p);
                    T_out(l,i)  = 1/W*Alpha(l,i,rho+1)*Beta(l,i)*d_V(l,o)*d_H(l,o)/s^2;
                end
            end
        end
    end
    T = T_FM(:,:,rho+1) + T_W(:,:,rho+1) + T_SP + T_SA + T_out;
    M_T = M_FM+M_SP+M_SA+M_out+M_W_SA;
    M_del = M_BRAM-M_T;
    [M_min, M_ind] = min(M_del);
    M_ind = M_ind+L*(0:(O*P*Q-1));
    T_at_min = T(M_ind);
    n_dsp_sel = 0;
    n_dsp_sel1 = 0;
    for n=1:O*P*Q
        if(n_dsp_sel<N_DSP)&&(n_dsp(n)<N_DSP)
            n_dsp_sel2 = n_dsp_sel1;
            n_dsp_sel1=n_dsp_sel;
            n_dsp_sel=n_dsp(n);
        end
    end
    DSP_selected2 = (n_dsp==n_dsp_sel2);
    DSP_selected1 = (n_dsp==n_dsp_sel1);
    DSP_selected = (n_dsp==n_dsp_sel);
    M_R = M_min.*(M_min>0).*DSP_selected;
    T_R = T_at_min.*(M_min>0).*DSP_selected;
    M_R1 = M_min.*(M_min>0).*DSP_selected1;
    T_R1 = T_at_min.*(M_min>0).*DSP_selected1;
    M_R2 = M_min.*(M_min>0).*DSP_selected2;
    T_R2 = T_at_min.*(M_min>0).*DSP_selected2;
    [T_R_v, T_R_i] = max((T_R==min(T_R(T_R>0)))*min(T_R(T_R>0)));
    [T_R_v1, T_R_i1] = max((T_R1==min(T_R1(T_R1>0)))*min(T_R1(T_R1>0)));
    [T_R_v2, T_R_i2] = max((T_R2==min(T_R2(T_R2>0)))*min(T_R2(T_R2>0)));
    if(T_R_i1<T_R_i)
        T_R_i = T_R_i1;
        T_R_v = T_R_v1;
        n_dsp_sel = n_dsp_sel1;
        M_R = M_R1;
    elseif(T_R_i2<T_R_i)
        T_R_i = T_R_i2;
        T_R_v = T_R_v2;
        n_dsp_sel = n_dsp_sel2;
        M_R = M_R2;
    end
    subplot(2,4,4*rho+1)
    bar(1:L,M_T(M_ind(T_R_i):M_ind(T_R_i)+L-1))
    xlabel('Layer Number');ylabel('\it M_{T}(i,l,\rho)(Mb)');
    title(strcat('c_{sa} = ',num2str(A(T_R_i)),', r_{sa} = ',num2str(B(T_R_i)),', n_{dsp} = ',num2str(n_dsp(T_R_i))))
    subplot(2,4,4*rho+2),scatter(n_dsp,M_min,Marker_Size,'p','filled','LineWidth', Line_Width);
    hold on,scatter(n_dsp_sel,M_min(T_R_i),5*Marker_Size,'p','filled','LineWidth', 3*Line_Width);
    hold on,set(gca,'xscale','log','FontSize', Font_Size),xticks(unique(n_dsp,'sorted'))
    txt = strcat('(',num2str(n_dsp_sel),',',num2str(M_min(T_R_i)));text(n_dsp_sel,M_min(T_R_i),txt)
    xlabel('DSP Units used');ylabel('\it M_{\Delta}(i,l,\rho)(Mb)');
    subplot(2,4,4*rho+3),scatter(n_dsp,Mb*T_at_min,Marker_Size,'filled','p','LineWidth', Line_Width);
    hold on,scatter(n_dsp_sel,Mb*T_R_v,5*Marker_Size,'filled','p','LineWidth', 3*Line_Width);
    hold on,set(gca,'xscale','log','FontSize', Font_Size),xticks(unique(n_dsp,'sorted'))
    txt = strcat('(',num2str(n_dsp_sel),',',num2str(Mb*T_R_v));text(n_dsp_sel,Mb*T_R_v,txt)
    xlabel('DSP Units used');ylabel('\it T(i) (million cycles)');
    subplot(2,4,4*rho+4),scatter3(n_dsp,M_min,Mb*T_at_min,Marker_Size,'p','filled','LineWidth', Line_Width);
    hold on,scatter3(n_dsp_sel,M_R(T_R_i),Mb*T_R_v,5*Marker_Size,'p','filled','LineWidth', 3*Line_Width);
    hold on,set(gca,'xscale','log','FontSize', Font_Size),xticks(unique(n_dsp,'sorted'))
    xlabel('DSP Units used');zlabel('\it T(i) (million cycles)');ylabel('\it M_{\Delta}(i,l,\rho)(Mb)');
    txt = strcat('(',num2str(n_dsp_sel),',',num2str(M_R(T_R_i)),',',num2str(Mb*T_R_v));text(n_dsp_sel,M_R(T_R_i),Mb*T_R_v,txt)
end