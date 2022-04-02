%MRC与CF的对比
%再加上PARA和SUCC

%用户和AP数目的规定
L=10;
M=100;

%循环次数
h_value=100;


%存放sum_rate
sum_rate_cf=zeros(h_value,1);
sum_rate_mrc_1=zeros(h_value,1);
sum_rate_mmse=zeros(h_value,1);

P_total=0.2*L;

%对不同的信道环境进行循环
for count_h=1:h_value
    
%信道矩阵的设置
%M个AP和K个用户，随机分布
%用户的数目
K=L; 
%分布的区域的大小
D= 0.25 ;
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
%带宽，Mhz
B=20;
%噪声功率
noise_p=B*10^6*1.381*10^(-23)*290*10^0.9;

%CF
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
end
[L1,M1]=size(A);
%标记A矩阵中元素为0的地方
flag_matrix=ones(M1,L1);
for m=1:M1
    for l=1:L1
        if A(l,m)==0
            flag_matrix(m,l)=0;
        end
    end
end
temp_rate1=zeros(L1,1);
for m=1:M1
a=A(:,m);
h=H_after(m,:);
local_0=find(a==0);
L1m=L1-length(local_0);
a(local_0)=[];
h(local_0')=[];
%计算temp_rate_1
M_temp=eye(L1m)-powervec_cf*h'*h/(powervec_cf*power(norm(h'),2)+1);
temp_rate1(m,1)=1/real(a'*M_temp*a);
end
[temp_rate1,temp_rate1_ind]=sort(temp_rate1,'descend');
rate_cf=zeros(1,L1);
for l=1:L1
    loca_1=find(flag_matrix(temp_rate1_ind,l)==1);
    loca_1=temp_rate1_ind(loca_1);
    for l1=1:M1
        if sum(ismember(temp_rate1_ind(l1),loca_1))~=0
           rate_cf(1,l)=max(0,log2(temp_rate1(l1)));
           break;
        end
    end
end
sum_rate_cf(count_h,1)=sum(rate_cf);


%MRC
R_mrc=zeros(L1,1);
for l=1:L1
    sum_mrc=0;
    for l1=1:L1
        h=H_after(:,l1)';
        if l1~=l
        sum_mrc=sum_mrc+norm(h*h')^2;
        end
    end
    R_mrc(l,1)=log2(1+(powervec_cf*norm(H_after(:,l)'*H_after(:,l))^2)/(norm(H_after(:,l))^2+powervec_cf*sum_mrc));
end      
sum_rate_mrc_1(count_h,1)=sum(R_mrc);

%MMSE
power_l=0.2;
noise_p_matrix=noise_p*eye(M);
rate_mmse=zeros(L1,1);
for l=1:L1
    sum_mmse=0;
    for l1=1:L1
        if l1~=l
            sum_mmse=sum_mmse+power_l*H_orig(:,l1)*H_orig(:,l1)';
        end
    end
    temp_mmse=abs(log2(1+(power_l*H_orig(:,l)'*inv(sum_mmse+noise_p_matrix)*H_orig(:,l))));
    rate_mmse(l)=temp_mmse;
end 
sum_rate_mmse(count_h,1)=sum(rate_mmse); 
count_h
end
mean(sum_rate_cf)
mean(sum_rate_mrc_1)
mean(sum_rate_mmse)