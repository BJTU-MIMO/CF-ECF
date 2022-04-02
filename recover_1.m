%tu 1
%首先检查候选集的有效性
%用户和AP数目的规定
L=10;
M_set=20:20:100;

%循环次数
h_value=30;
M_value=5;


%存放sum_rate
sum_rate_cf_part=zeros(h_value,M_value);
sum_rate_cf_all=zeros(h_value,M_value);
sum_rate_para_all=zeros(h_value,M_value);
sum_rate_para_part=zeros(h_value,M_value);
sum_rate_mmse=zeros(h_value,M_value);

P_total=0.2*L;

for count_M=1:M_value
M=M_set(count_M);
% M=100;
%对不同的信道环境进行循环
for count_h=1:h_value

[H_orig,BETAA] = channel(L,M);
%带宽，Mhz
B=20;
%噪声功率
noise_p=B*10^6*1.381*10^(-23)*290*10^0.9;
powervec=((P_total/L)/noise_p)*ones(1,L);
L1=L;
M1=M;
[L1,M1,A,H_after,flag_matrix] = A_matrix(L1,M1,powervec,H_orig );
A_cf=A;
H_cf=H_after;
if rank(A)<L1
    continue;
end


%cf
temp_rate1=zeros(M1,1);
sigema_cf_1=zeros(M1,1);
powervec_cf=((P_total/L)/noise_p);
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
sigema_cf_1(m,1)=powervec_cf/temp_rate1(m,1);
end
[temp_rate1,temp_rate1_ind_1]=sort(temp_rate1,'descend');
rate_cf_1=zeros(1,L1);
%记录到底哪些AP被选择
AP_selec=zeros(M1,L1);
for l=1:L1
    loca_1=find(flag_matrix(temp_rate1_ind_1,l)==1);
    loca_1=temp_rate1_ind_1(loca_1);
    for l1=1:M1
        if sum(ismember(temp_rate1_ind_1(l1),loca_1))~=0
           rate_cf_1(1,l)=max(0,log2(temp_rate1(l1)));
           AP_selec(temp_rate1_ind_1(l1),l)=1;
           break;
        end
    end
end
[~,ind_AP]=find(sum(transpose(AP_selec))~=0);
sum_rate_cf_all(count_h,count_M)=sum(rate_cf_1);

%挑选出候选集合
[ flag_matrix_hm,A_hm_final,H_hm_final ] = APselect( flag_matrix,H_after,A,L1,M1);

%cf_hm
temp_rate1=zeros(L1,1);
powervec_cf=((P_total/L)/noise_p);
for m=1:L1
a=A_hm_final(:,m);
h=H_hm_final(m,:);
local_0=find(a==0);
L1m=L1-length(local_0);
a(local_0)=[];
h(local_0')=[];
%计算temp_rate_1
M_temp=eye(L1m)-powervec_cf*h'*h/(powervec_cf*power(norm(h'),2)+1);
temp_rate1(m,1)=1/real(a'*M_temp*a);
end
[temp_rate1,temp_rate1_ind_1]=sort(temp_rate1,'descend');
rate_cf_2=zeros(1,L1);
for l=1:L1
    loca_1=find(flag_matrix_hm(temp_rate1_ind_1,l)==1);
    loca_1=temp_rate1_ind_1(loca_1);
    for l1=1:L1
        if sum(ismember(temp_rate1_ind_1(l1),loca_1))~=0
           rate_cf_2(1,l)=max(0,log2(temp_rate1(l1)));
           break;
        end
    end
end
sum_rate_cf_part(count_h,count_M)=sum(rate_cf_2);

%all_para
%SOCP规划
%每一次迭代之后还要对a进行迭代
flag_para_all=0;
t_max=powervec_cf*L1;
t_min=10^9;
epsi=10^4;
A_para=A;
H_para=H_after;
while(t_max-t_min>epsi)
A_abs=zeros(L1,M1);
H_abs=zeros(M1,L1);
for l=1:L1
    for m=1:M1
        A_abs(l,m)=abs(A_para(l,m));
        H_abs(m,l)=abs(H_para(m,l));
    end
end
% Q=zeros(L1,L1,M1);
% for m=1:M1
%     for l1=1:L1
%         for l2=1:L1
%             Q(l1,l2,m)=2*A_abs(l1,m)*A_abs(l2,m)*H_abs(m,l1)*H_abs(m,l2);
%         end
%     end
% end          
t=(t_max+t_min)/2;
powervec_para_sdpvar=sdpvar(L1,1);
F=[
    sum(powervec_para_sdpvar)==powervec_cf*L1
    ];
for l=1:L1
    F=F+[
        powervec_para_sdpvar(l,1)>=0
        ];
end
for m=1:M1
%     F=F+[
%         1+sum(powervec_para_sdpvar'.*H_abs(m,:).*H_abs(m,:))<=(1/t)*(1/2)*powervec_para_sdpvar'*Q(:,:,m)*powervec_para_sdpvar
%         ];
     F=F+[
         norm(sqrt(powervec_para_sdpvar).* H_abs(m,:)')<=(1/sqrt(t))*powervec_para_sdpvar.*H_abs(m,:)'.*A_abs(:,m)
         ];
end
% for m=1:M1
%     F=F+[
%         sum(powervec_para_sdpvar.*A_abs(:,m).*A_abs(:,m))-t<=powervec_cf*L1       
%         ];
% end
ops=sdpsettings('solver','mosek+','verbose',1,'cachesolvers',0);
result=optimize(F,0,ops)
if result.problem==0
    t_min=t;
    powervec_para_all=double(powervec_para_sdpvar);
    [L1,M1,A_para,H_after,flag_matrix] = A_matrix(L1,M1,powervec_para_all,H_orig );
    H_para=H_after;
    flag_para_all=flag_para_all+1;
else
    t_max=t;  
end
if t_max-t_min<epsi
    break;
end
end
if flag_para_all==0
    continue;
end
%计算rate
% 按照每个用户真正参与的线性组合计算有效噪声以及SE
sigema_para_all=zeros(M1,1);
for m=1:M1
a=A_para(:,m);
h=H_para(m,:);
local_0=find(a==0);
L1m=L1-length(local_0);
a(local_0)=[];
h(local_0')=[];
p=powervec_para_all;
p(local_0')=[];
% 计算有效噪声
sigema_para_all(m)=real(a'*(diag(p)-(1/(1+h*diag(p)*h'))*diag(p)*h'*h*diag(p))*a);
end
% 每个用户找到自己参与的线性组合，取最小的噪声，还要保证A矩阵满秩
rate_para=zeros(1,L1);
[sigema_para_all,sigema_para_ind]=sort(sigema_para_all);
sigema_local_0=find(sigema_para_all==0);
sigema_para_all(sigema_local_0)=[];
sigema_length=length(sigema_para_all);
sigema_para_ind(sigema_local_0)=[];
for l=1:L1
    loca_1=find(flag_matrix(sigema_para_ind,l)==0);
    loca_1=sigema_para_ind(loca_1);
    for l1=1:sigema_length
        if sum(ismember(sigema_para_ind(l1),loca_1))~=0
            rate_para(1,l)=max(0,log2(abs(powervec_para_all(l)/sigema_para_all(l1))));
            break;
        end
    end
end
sum_rate_para_all(count_h,count_M)=sum(rate_para);


%part_para
%二次规划
%每一次迭代之后还要对a进行迭代
flag_para_part=0;
t_max=powervec_cf*L1;
t_min=10^9;
epsi=10^4;
A_para_part=A_hm_final;
H_para_part=H_hm_final;
while(t_max-t_min>epsi)
A_abs=zeros(L1,L1);
H_abs=zeros(L1,L1);
for l=1:L1
    for m=1:L1
        A_abs(l,m)=abs(A_para_part(l,m));
        H_abs(m,l)=abs(H_para_part(m,l));
    end
end
Q=zeros(L1,L1,L1);
for m=1:L1
    for l1=1:L1
        for l2=1:L1
            Q(l1,l2,m)=2*A_abs(l1,m)*A_abs(l2,m)*H_abs(m,l1)*H_abs(m,l2);
        end
    end
end          
t=(t_max+t_min)/2;
powervec_para_sdpvar=sdpvar(L1,1);
F=[
    sum(powervec_para_sdpvar)==powervec_cf*L1
    ];
for l=1:L1
    F=F+[
        powervec_para_sdpvar(l,1)>=0
        ];
end
for m=1:L1
    F=F+[
        (1/2)*powervec_para_sdpvar'*Q(:,:,m)*powervec_para_sdpvar-sum(powervec_para_sdpvar'.*H_abs(m,:).*H_abs(m,:))*t>=t
        ];
end
for m=1:L1
    F=F+[
        sum(powervec_para_sdpvar.*A_abs(:,m).*A_abs(:,m))-t<=powervec_cf*L1       
        ];
end
ops=sdpsettings('solver','mosek+','verbose',1,'cachesolvers',0);
result=optimize(F,0,ops)
if result.problem==0
    t_min=t;
    powervec_para_part=double(powervec_para_sdpvar);
    [L1,L1,A_para_part,H_para_part,flag_matrix_hm] = A_matrix(L1,L1,powervec_para_part,H_para_part );
    flag_para_part=flag_para_part+1;
else
    t_max=t;  
end
if t_max-t_min<epsi
    break;
end
end
if flag_para_part==0
    continue;
end
%计算rate
% 按照每个用户真正参与的线性组合计算有效噪声以及SE
sigema_para_part=zeros(L1,1);
for m=1:L1
a=A_para_part(:,m);
h=H_para_part(m,:);
local_0=find(a==0);
L1m=L1-length(local_0);
a(local_0)=[];
h(local_0')=[];
p=powervec_para_part;
p(local_0')=[];
% 计算有效噪声
sigema_para_part(m)=real(a'*(diag(p)-(1/(1+h*diag(p)*h'))*diag(p)*h'*h*diag(p))*a);
end
% 每个用户找到自己参与的线性组合，取最小的噪声，还要保证A矩阵满秩
rate_para_part=zeros(1,L1);
[sigema_para_part,sigema_para_ind]=sort(sigema_para_part);
sigema_local_0=find(sigema_para_part==0);
sigema_para_part(sigema_local_0)=[];
sigema_length=length(sigema_para_part);
sigema_para_ind(sigema_local_0)=[];
for l=1:L1
    loca_1=find(flag_matrix_hm(sigema_para_ind,l)==0);
    loca_1=sigema_para_ind(loca_1);
    for l1=1:sigema_length
        if sum(ismember(sigema_para_ind(l1),loca_1))~=0
            rate_para_part(1,l)=max(0,log2(abs(powervec_para_part(l)/sigema_para_part(l1))));
            break;
        end
    end
end
sum_rate_para_part(count_h,count_M)=sum(rate_para_part);


%MMSE
power_l=P_total/L;
noise_p_matrix=noise_p*eye(M);
rate_mmse=zeros(L,1);
for l=1:L
    sum_mmse=0;
    for l1=1:L
        if l1~=l
            sum_mmse=sum_mmse+power_l*H_orig(:,l1)*H_orig(:,l1)';
        end
    end
    temp_mmse=abs(log2(1+(power_l*H_orig(:,l)'*pinv(sum_mmse+noise_p_matrix)*H_orig(:,l))));
    rate_mmse(l)=temp_mmse;
end 
sum_rate_mmse(count_h,count_M)=sum(rate_mmse);


count_h
end
count_M
end