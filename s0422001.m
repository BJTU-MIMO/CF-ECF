%Fig.1
%修改用户数和AP数的比例
L=40;
M=100;


%循环次数
h_value=100;



%存放sum_rate
sum_rate_cf=zeros(h_value,1);
sum_rate_para=zeros(h_value,1);
sum_rate_hm_para=zeros(h_value,1);


P_total=0.1*L;
%对不同的信道环境进行循环
for count_h=1:h_value
    
%信道矩阵的设置
%M个AP和K个用户，随机分布
%用户的数目
K=L; 
%分布的区域的大小
D=0.2;
%基站天线的高度
Hb = 15; 
%移动台天线的高度
Hm = 1.65; 
%载频，MHz
f = 1900; 
aL = (1.1*log10(f)-0.7)*Hm-(1.56*log10(f)-0.8);
L_loss = 46.3+33.9*log10(f)-13.82*log10(Hb)-aL;
%参考距离
d0=0.01;%km
d1=0.05;%km
%Randomly locations of M APs:
AP=unifrnd(-D/2,D/2,M,2);
%Wrapped around (8 neighbor cells)
D1=zeros(M,2);
D1(:,1)=D1(:,1)+ D*ones(M,1);
AP1=AP+D1;
D2=zeros(M,2);
D2(:,2)=D2(:,2)+ D*ones(M,1);
AP2=AP+D2;
D3=zeros(M,2);
D3(:,1)=D3(:,1)- D*ones(M,1);
AP3=AP+D3;
D4=zeros(M,2);
D4(:,2)=D4(:,2)- D*ones(M,1);
AP4=AP+D4;
D5=zeros(M,2);
D5(:,1)=D5(:,1)+ D*ones(M,1);
D5(:,2)=D5(:,2)- D*ones(M,1);
AP5=AP+D5;
D6=zeros(M,2);
D6(:,1)=D6(:,1)- D*ones(M,1);
D6(:,2)=D6(:,2)+ D*ones(M,1);
AP6=AP+D6;
D7=zeros(M,2);
D7=D7+ D*ones(M,2);
AP7=AP+D7;
D8=zeros(M,2);
D8=D8- D*ones(M,2);
AP8=AP+D8;
%Randomly locations of K terminals:
Ter=unifrnd(-D/2,D/2,K,2);
sigma_shd=8; %in dB
%Create an MxK large-scale coefficients beta_mk
BETAA = zeros(M,K);
betadB=zeros(M,K);
dist=zeros(M,K);
%路径损耗和小尺度衰落
h_s=sqrt(1/2)*(randn(M,K) +1j*randn(M,K));
%h计算并存放信道系数
H_orig = zeros(M,K);
for m=1:M
    for k=1:K
    dist(m,k) = min([norm(AP(m,:)-Ter(k,:)), norm(AP1(m,:)-Ter(k,:)),norm(AP2(m,:)-Ter(k,:)),norm(AP3(m,:)-Ter(k,:)),norm(AP4(m,:)-Ter(k,:)),norm(AP5(m,:)-Ter(k,:)),norm(AP6(m,:)-Ter(k,:)),norm(AP7(m,:)-Ter(k,:)),norm(AP8(m,:)-Ter(k,:)) ]); %distance between Terminal k and AP m
    if dist(m,k)<d0
         betadB(m,k)=-L_loss - 35*log10(d1) + 20*log10(d1) - 20*log10(d0);
    elseif ((dist(m,k)>=d0) && (dist(m,k)<=d1))
         betadB(m,k)= -L_loss - 35*log10(d1) + 20*log10(d1) - 20*log10(dist(m,k));
    else
    betadB(m,k) = -L_loss - 35*log10(dist(m,k)) + sigma_shd*randn(1,1); %large-scale in dB
    end
    BETAA(m,k)=10^(betadB(m,k)/10); 
    H_orig(m,k)=sqrt(BETAA(m,k))*h_s(m,k);
    end
end
BETAA_g=BETAA;
%带宽，Mhz
B=20;
%噪声功率
noise_p=B*10^6*1.381*10^(-23)*290*10^0.9;

%求系数矩阵A
%存放系数A
A=zeros(L,M);
A_real=zeros(L,M);
A_imag=zeros(L,M);
%功率平分
powervec_cf=(P_total/L)/noise_p;
%对不同的AP进行循环
for m=1:M
%信道矩阵分成实部和虚部两部分
%g_m以列向量的形式参与运算
H_temp=real(H_orig(m,:))';
for count_part=1:2
%第二步
%计算ak_
[h,perm_part]=sort(abs(H_temp));
temp_k=sqrt(1 + powervec_cf*(h'*h));
K=floor(temp_k);
%计算a1_
G=eye(L)-powervec_cf*(h*h')/(1 + powervec_cf*(h'*h));
r=-pinv(G(1:L-1,1:L-1))*G(1:L-1,L);
a1_=[r;1];
%计算ak_
a=zeros(L,K);
ak_=mat2cell(a,L,ones(1,K));
temp_ak_=zeros(L,1);
for k=1:K
    temp_ak_=cell2mat(ak_(1,k));
    for l=1:L
        temp_ak_(l,1)=k*a1_(l,1);
    end
    ak_(1,k)={temp_ak_}; 
end
%第三步，successive quantization
%按照自己的理解修改
temp2=zeros(L,1);
for k=1:K
    for l=1:L
        temp1=cell2mat(ak_(1,k));
        if l==L
            temp2=temp1;
            break;
        end
        temp1(l,1)=floor(temp1(l,1));
        temp3=temp1;
        temp1=cell2mat(ak_(1,k));
        temp1(l,1)=ceil(temp1(l,1));
        temp4=temp1;
        temp5=transpose(temp3)*G*temp3;
        temp6=transpose(temp4)*G*temp4;
        if temp5<temp6
            temp2=temp3;
        else
            temp2=temp4;
        end
        ak_(1,k)={temp2};
    end
end
%第四步
temp_a_final1=zeros(K,1);
a_final=zeros(L,1);
for k=1:K
    temp_a_final1(k,1)=transpose(cell2mat(ak_(1,k)))*G*cell2mat(ak_(1,k));
end
[value,temp_a_final2]=min(temp_a_final1);
for l=1:L
    temp_a_final3=cell2mat(ak_(1,temp_a_final2));
    a_final(l,1)=temp_a_final3(l,1);
end
A_final=zeros(L,1);
l_A_count_1=0;
for l_A_count_2=1:L
    while 1
        l_A_count_1=l_A_count_1+1;
        if l_A_count_1==perm_part(l_A_count_2,1)
            break;
        end
    end
    A_final(l_A_count_2,1)=a_final(perm_part(l_A_count_2,1),1);
    l_A_count_1=0;
end
if count_part==1
    A_real(:,m)=A_final;
else
    A_imag(:,m)=A_final;
end       
A(:,m)=A_real(:,m)+1j*A_imag(:,m);
H_temp=imag(H_orig(m,:))';
end
H_after=abs(real(H_orig))+1j*abs(imag(H_orig));
end
%判断A中有无为a_m=0
ind_m=zeros(1,M);
for m=1:M
   if A(:,m)==zeros(L,1)
      ind_m(1,m)=1; 
   end
end
if ind_m~=zeros(1,M)
A(:,ind_m)=[];
H_after(ind_m',:)=[];
BETAA_g(ind_m',:)=[];
end
[L1,M1]=size(A);
%判断A中有无为a_l=0
ind_l=zeros(L1,1);
for l=1:L1
   if A(l,:)==zeros(1,M1)
      ind_l(l)=1; 
   end
end
if ind_l~=zeros(L1,1)
A(ind_l,:)=[];
H_after(:,ind_l')=[];
BETAA_g(:,ind_l')=[];
end
[L1,M1]=size(A);
%标记A矩阵中元素为0的地方
flag_matrix=zeros(M1,L1);
for m=1:M1
    for l=1:L1
        if A(l,m)==0
            A(l,m)=10^(-4)+1j*10^(-4);
            flag_matrix(m,l)=1;
        end
    end
end


%计算几何均值
BETAA_mean=zeros(M1,1);
for m=1:M1
    geo_mean_temp=abs(A(:,m)').*BETAA_g(m,:);
    BETAA_mean(m,1)=geo_mean(geo_mean_temp);
end
%排序，挑选候选的子集
[~,BETAA_permu_ind]=sort(BETAA_mean,'ascend');
A_para=A(:,BETAA_permu_ind');
H_para=H_after(BETAA_permu_ind,:);
flag_matrix=flag_matrix(BETAA_permu_ind,:);
A_para_final=A_para(:,1);
ind_a=zeros(L1,1);
ind_a(1)=1;
count_ind=1;
%挑出满秩的部分
for m=2:M1
    a=[A_para_final,A_para(:,m)];
    if rank(a')==size(a,2)
        A_para_final=a;
        count_ind=count_ind+1;
        ind_a(count_ind)=m;
    end
    if rank(a')==L1
        break;
    end
end
H_para_final=H_after(ind_a,:);
H_para_final_a=abs(H_para_final);
A_para_final_a=abs(A_para_final);
%判断用户到底被哪几个AP服务
flag_matrix_para=flag_matrix(ind_a,:);

%普通的cf,没有候选集合
temp_rate1=zeros(L1,1);
%存放sigema_cf
sigema_cf=zeros(M1,1);
for m=1:M1
a=A(:,m);
h=H_after(m,:);
%计算temp_rate_1
M_temp=eye(L1)-powervec_cf*h'*h/(powervec_cf*power(norm(h'),2)+1);
temp_rate1(m,1)=real(a'*M_temp*a);
sigema_cf(m,1)=powervec_cf*temp_rate1(m,1);
end
%每个用户找到自己参与的线性组合，取最小的噪声，还要保证A矩阵满秩
rate_cf=zeros(1,L1);
%升序排序
[sigema_cf,sigema_cf_ind]=sort(sigema_cf);
A_cf_final=A(:,sigema_cf_ind');
H_cf_final=H_after(sigema_cf_ind,:);
for l=1:L1
    loca_1=find(flag_matrix(sigema_cf_ind,l)==0);
    loca_1=sigema_cf_ind(loca_1);
    for l1=1:M1
        if sum(ismember(sigema_cf_ind(l1),loca_1))~=0
           rate_cf(1,l)=max(0,log2(abs(powervec_cf/sigema_cf(l1))));
           break;
        end
    end
end
sum_rate_cf(count_h,1)=sum(rate_cf);

% para
t_max=10^12;
t_min=10^9;
epsi=10^3;
while(t_max-t_min>epsi)
t=(t_max+t_min)/2;
powervec_para_sdpvar=sdpvar(1,L1);
r_m=sdpvar(1,1);
s_m=sdpvar(1,1);
F=[
    sum(powervec_para_sdpvar)==powervec_cf*L1
    ];
for m=1:L1
    F=F+[
        powervec_para_sdpvar(1,L1)>=0
        ];
end
for m1=1:L1
    F=F+[
    sum((((A_para_final_a(:,m)').^2)-conj(H_para_final_a(m,:)).*conj(A_para_final_a(:,m)')*r_m+((A_para_final_a(:,m)').^2)*s_m).*powervec_para_sdpvar)<=t,
    sum(powervec_para_sdpvar.*H_para_final_a(m,:).*A_para_final_a(:,m)')<=r_m,
    sum(powervec_para_sdpvar.*(H_para_final_a(m,:).^2))<=s_m
    ];
end
ops=sdpsettings('solver','','verbose',0,'cachesolvers',0);
result=optimize(F,0,ops);
if result.problem==0
    t_max=t;
    powervec_para=double(powervec_para_sdpvar);
else
     t_min=t;
end
if t_max-t_min<epsi
    break;
end
end
% 按照每个用户真正参与的线性组合计算有效噪声以及SE
sigema_para=zeros(L1,1);
for m=1:L1
a=A_para_final(:,m);
h=H_para_final(m,:);
powervec_para_m=powervec_para;
% 计算有效噪声
sigema_para(m)=real(a'*(diag(powervec_para_m)-(1/(1+h*diag(powervec_para_m)*h'))*diag(powervec_para_m)*h'*h*diag(powervec_para_m))*a);
end
% 每个用户找到自己参与的线性组合，取最大的噪声，还要保证A矩阵满秩
rate_para=zeros(1,L1);
[sigema_para,sigema_para_ind]=sort(sigema_para);
for l=1:L1
    loca_1=find(flag_matrix_para(sigema_para_ind,l)==0);
    loca_1=sigema_para_ind(loca_1);
    for l1=1:L1
        if sum(ismember(sigema_para_ind(l1),loca_1))~=0
            rate_para(1,l)=max(0,log2(abs(powervec_para(l)/sigema_para(l1))));
            break;
        end
    end
end
sum_rate_hm_para(count_h,1)=sum(rate_para);

%无候选集合
% para
t_max=10^12;
t_min=10^9;
epsi=10^3;
while(t_max-t_min>epsi)
t=(t_max+t_min)/2;
powervec_para_sdpvar=sdpvar(1,L1);
r_m=sdpvar(1,1);
s_m=sdpvar(1,1);
F=[
    sum(powervec_para_sdpvar)==powervec_cf*L1
    ];
for m=1:L1
    F=F+[
        powervec_para_sdpvar(1,L1)>=0
        ];
end
for m1=1:M1
    F=F+[
    sum((((A(:,m)').^2)-conj(H_after(m,:)).*conj(A(:,m)')*r_m+((A(:,m)').^2)*s_m).*powervec_para_sdpvar)<=t,
    sum(powervec_para_sdpvar.*H_after(m,:).*A(:,m)')<=r_m,
    sum(powervec_para_sdpvar.*(H_after(m,:).^2))<=s_m
    ];
end
ops=sdpsettings('solver','','verbose',0,'cachesolvers',0);
result=optimize(F,0,ops);
if result.problem==0
    t_max=t;
    powervec_para=double(powervec_para_sdpvar);
else
     t_min=t;
end
if t_max-t_min<epsi
    break;
end
end
% 得到P之后的调整
% 按照每个用户真正参与的线性组合计算有效噪声以及SE
sigema_para=zeros(M1,1);
for m=1:M1
a=A(:,m);
h=H_after(m,:);
powervec_para_m=powervec_para;
% 计算有效噪声
sigema_para(m)=real(a'*(diag(powervec_para_m)-(1/(1+h*diag(powervec_para_m)*h'))*diag(powervec_para_m)*h'*h*diag(powervec_para_m))*a);
end
% 每个用户找到自己参与的线性组合，取最大的噪声，还要保证A矩阵满秩
rate_para=zeros(1,L1);
[sigema_para,sigema_para_ind]=sort(sigema_para);
for l=1:L1
    loca_1=find(flag_matrix(sigema_para_ind,l)==0);
    loca_1=sigema_para_ind(loca_1);
    for l1=1:M1
        if sum(ismember(sigema_para_ind(l1),loca_1))~=0
            rate_para(1,l)=max(0,log2(abs(powervec_para(l)/sigema_para(l1))));
            break;
        end
    end
end
sum_rate_para(count_h,1)=sum(rate_para);

count_h
end


