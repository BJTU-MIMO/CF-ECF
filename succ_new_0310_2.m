%���г�����д
%ǰ���ŵ������Լ�����cf��A����Ĳ�����ȷ
%�Լ���Aʱ��quantization���Լ�������޸�
%cf�޸����
%para���Ż�������,��ʱ��Ϊ��������쿴���ѭ���Ľ��
%������cf,para,succ��ֻ��ҪL���������
%δ�Ż���succ����д���ˣ�δ����
%��ȷ��para�Ż�������ȷ�ԣ��Ż�����succ,��ȷ�Ի�������ȷ��������Ϊʲô���ж���ô��Σ�
%succ��beta�㷨�����⣬�Ѿ��޸ģ�ȷ���û�����˳���rate��Ӱ��
%����д��
%����������
%ȥ���Ż���beta���㷨�����Ż�ǰ����
%betaд�������⣬��L=2ʱ��succ_2Ӧ��3��4�е�һ�����,�ѽ��
%�Ż�beta�㷨


%�޸��û�����AP���ı���
L=2;
M=4;

%����ȡ���ٸ�ֵ
x_value=3;

%ѭ������
h_value=200;

%���sum_rate
sum_rate_cf=zeros(h_value,x_value);
sum_rate_para=zeros(h_value,x_value);
%�ȳ����ڼ���succ��ʱ�������para�Ż�����P����
%��succ�����Ż��ĳ��������
%����succ���й����Ż�
%succ_1 ����ȫ���еķ���ȷ���û�����˳��
sum_rate_succ_1=zeros(h_value,x_value);
%succ_3 �û�������˳��ͨ��betaȷ����betaֵ����û��Ƚ������룬�Ա���succ_1��������
sum_rate_succ_3=zeros(h_value,x_value);
%succ_4 �û�������˳��ͨ��betaȷ����betaֵС���û��Ƚ������룬�Ա���succ_2�������٣�ȷ���û�������˳������Ӱ��
sum_rate_succ_4=zeros(h_value,x_value);
%��succ�����Ż�
%succ_2 �û�������˳��ͨ��ȫ����Ѱ��
sum_rate_succ_2=zeros(h_value,x_value);


%�趨final����ȥ��rateΪ0�����
final_sum_rate_cf=zeros(1,x_value);
final_sum_rate_para=zeros(1,x_value);
final_sum_rate_succ_1=zeros(1,x_value);
final_sum_rate_succ_2=zeros(1,x_value);
final_sum_rate_succ_3=zeros(1,x_value);
final_sum_rate_succ_4=zeros(1,x_value);

%���per_rate
per_rate_cf=zeros(h_value,x_value);
per_rate_para=zeros(h_value,x_value);
per_rate_succ_1=zeros(h_value,x_value);
per_rate_succ_2=zeros(h_value,x_value);
per_rate_succ_3=zeros(h_value,x_value);
per_rate_succ_4=zeros(h_value,x_value);


final_per_rate_cf=zeros(1,x_value);
final_per_rate_para=zeros(1,x_value);
final_per_rate_succ_1=zeros(1,x_value);
final_per_rate_succ_2=zeros(1,x_value);
final_per_rate_succ_3=zeros(1,x_value);
final_per_rate_succ_4=zeros(1,x_value);


%���ʵĳ�ʼֵ���ܹ���
%���ܹ��ʽ��б仯
P_total=zeros(1,x_value);


%�Բ�ͬ���ŵ���������ѭ��
for count_h=1:h_value
 
%�ŵ����������
%M��AP��K���û�������ֲ�
%�û�����Ŀ
K=L; 
%�ֲ�������Ĵ�С
D=1;
%��վ���ߵĸ߶�
Hb = 15; 
%�ƶ�̨���ߵĸ߶�
Hm = 1.65; 
%��Ƶ��MHz
f = 1900; 
aL = (1.1*log10(f)-0.7)*Hm-(1.56*log10(f)-0.8);
L_loss = 46.3+33.9*log10(f)-13.82*log10(Hb)-aL;
%�ο�����
d0=0.01;%km
d1=0.05;%km
%M��AP����ֲ�
%����һ��ά��ΪM*2�ľ��󣬾�����ÿ��Ԫ�ؾ�Ϊ(-D/2,D/2)����ȷֲ��������
AP=unifrnd(-D/2,D/2,M,2);
%Wrapped around (8 neighbor cells)
%��ֹ��ԵЧӦ�����а���(8�����ŵ�cell��
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
%K���û�����ֲ�
%����һ��ά��ΪK*2�ľ��󣬾�����ÿ��Ԫ�ص�ȡֵΪ(-D/2,D/2)����ȷֲ��������
Ter=unifrnd(-D/2,D/2,K,2);
%����·����ĵĲ�����in dB
sigma_shd=8; 
%·����ĺ�С�߶�˥��
z_s=randn(M,K);
p_s=10.^(sigma_shd*z_s/10);
h_s=sqrt(1/2)*(randn(M,K) +1j*randn(M,K));
%h���㲢����ŵ�ϵ��
H_orig = zeros(M,K);
%����һ����СΪM*K�Ĵ��beta�ľ���
BETAA = zeros(M,K);
%dist����û���AP��֮��ľ���
dist=zeros(M,K);
for m=1:M
    for k=1:K
    %�������������AP��ľ���
    dist(m,k) = min([norm(AP(m,:)-Ter(k,:)), norm(AP1(m,:)-Ter(k,:)),norm(AP2(m,:)-Ter(k,:)),norm(AP3(m,:)-Ter(k,:)),norm(AP4(m,:)-Ter(k,:)),norm(AP5(m,:)-Ter(k,:)),norm(AP6(m,:)-Ter(k,:)),norm(AP7(m,:)-Ter(k,:)),norm(AP8(m,:)-Ter(k,:)) ]); %distance between Terminal k and AP m
    if dist(m,k)<d0
         betadB=-L_loss - 35*log10(d1) + 20*log10(d1) - 20*log10(d0);
    elseif ((dist(m,k)>=d0) && (dist(m,k)<=d1))
         betadB= -L_loss - 35*log10(d1) + 20*log10(d1) - 20*log10(dist(m,k));
    else
         betadB = -L_loss - 35*log10(dist(m,k));
    end  
    BETAA(m,k)=10^(betadB/10); 
    H_orig(m,k)=sqrt(BETAA(m,k)*p_s(m,k))*h_s(m,k);
    end
end

for count_p=1:x_value
P_total(1,count_p)=(0.3+0.1*count_p)*L;
%����Mhz
B=20;
%��������
noise_p=B*10^6*1.381*10^(-23)*290*10^0.9;


%��ʼ����cf_rate
%���ϵ��A
A=zeros(L,M);
A_real=zeros(L,M);
A_imag=zeros(L,M);
%���sigema_cf
sigema_cf=zeros(M,1);
%��ʱ��Ų�ͬ����Ч�����������rate,ȡС
temp_rate1=zeros(M,1);
%����ƽ��
powervec_cf=(P_total(1,count_p)/L)/noise_p;

%�Բ�ͬ��AP����ѭ��
for m=1:M
%�ŵ�����ֳ�ʵ�����鲿������
%g_m������������ʽ��������
H_temp=real(H_orig(m,:))';
for count_part=1:2
%�ڶ���
%����ak_
[h,perm_part]=sort(abs(H_temp));
temp_k=sqrt(1 + powervec_cf*(h'*h));
K=floor(temp_k);
%����a1_
G=eye(L)-powervec_cf*(h*h')/(1 + powervec_cf*(h'*h));
r=-pinv(G(1:L-1,1:L-1))*G(1:L-1,L);
a1_=[r;1];
%����ak_
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
%��������successive quantization
%�����Լ�������޸�
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
%���Ĳ�
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
%����temp_rate_1
H_after=abs(real(H_orig))+1j*abs(imag(H_orig));
M_temp=eye(L)-powervec_cf*transpose(H_after(m,:))*conj(H_after(m,:))/(powervec_cf*power(norm(transpose(H_after(m,:))),2)+1);
temp_rate1(m,1)=real(ctranspose(A(:,m))*M_temp*A(:,m));
sigema_cf(m,1)=powervec_cf*temp_rate1(m,1);
end
[sigema_cf,sigema_cf_ind]=sort(sigema_cf);
rate_cf=max(0,log2(1/temp_rate1(sigema_cf_ind(L,1),1)));
per_rate_cf(count_h,count_p)=rate_cf;
sum_rate_cf(count_h,count_p)=per_rate_cf(count_h,count_p)*L;


%��ʼδ�������Ż���succ����
powervec_ecf=(P_total(1,count_p)/L)*ones(1,L)/noise_p;
P_succ_a=diag(powervec_ecf);
%��A������sigema�����˳���������
A_trans=transpose(A);
A_succ_a=A_trans(sigema_cf_ind,:);
flag1=rank(A_succ_a);
A_succ_a=A_succ_a(1:L,:);
%����H
H_succ_a=H_after(sigema_cf_ind,:);
H_succ_a=H_succ_a(1:L,:);

%,succ_rate_1 ȫ����
%ȫ���еķ���
%�ж��������з�ʽ
per_temp=perms(1:L);
[num1,~]=size(per_temp);
rate_succ_1=zeros(L,num1);
sigema_succ_1=zeros(L,num1);
for n=1:num1
    A_succ_1=A_succ_a(:,per_temp(n,:));
    H_succ_1=H_succ_a(:,per_temp(n,:));
    %��ʼ����rate
    for m=1:L
        if m==1
            sigema_succ_1(1,n)=A_succ_1(1,:)*pinv(pinv(P_succ_a)+transpose(H_succ_1(1,:))*conj(H_succ_1(1,:)))*A_succ_1(1,:)';
        else
            ins_chol_1=pinv(pinv(P_succ_a)+transpose(H_succ_1(m,:))*conj(H_succ_1(m,:)));
            ins_chol_2=ins_chol_1-ins_chol_1*A_succ_1(1:(m-1),:)'*pinv(A_succ_1(1:(m-1),:)')*pinv(ins_chol_1)*pinv(A_succ_1(1:(m-1),:))*A_succ_1(1:(m-1),:)*ins_chol_1;
            sigema_succ_1(m,n)=abs(real(A_succ_1(m,:)*ins_chol_2*A_succ_1(m,:)'));
        end
    end
    [sigema_succ_1_permu,~]=sort(sigema_succ_1(:,n));
    temp_L=per_temp(n,:);
    for count_temp_L=1:L
        rate_succ_1(temp_L(1,count_temp_L),n)=max(0,real(log2(powervec_ecf(1,count_temp_L)/sigema_succ_1_permu(count_temp_L,1))));
    end
end
sum_rate_succ_1_temp=sum(rate_succ_1);
[max_value,max_pos]=max(sum_rate_succ_1_temp);
per_rate_succ_1(count_h,count_p)=max(rate_succ_1(:,max_pos));
sum_rate_succ_1(count_h,count_p)=max_value;


%��ʼ����para_rate
%A���������ı�
%��ʼ��P�����Ż����Ӽ���sigema��ʼ��ȫ����д�ɷ��ź�������������֪�������ֱ�ӽ����Ż�����
syms P1 P2 value_1
powervec_para=[P1,P2];
P_para=diag(powervec_para);
sigema_para_h_temp=sym('H',[M,L]);
sigema_para_a_temp=sym('A',[L,M]);
rate_para_temp=sym('R',[M,L]);
sigema_para_temp=sym('sigema',[M,1]);
for m=1:M
    sigema_para_temp(m,1)=sigema_para_a_temp(:,m)'*pinv(pinv(P_para)+transpose(sigema_para_h_temp(m,:))*conj(sigema_para_h_temp(m,:)))*sigema_para_a_temp(:,m);
    for l=1:L
        rate_para_temp(m,l)=real(log2((P_para(l,l)/sigema_para_temp(m,1))));
    end
end
goal_para=value_1*sum(sum(rate_para_temp));
f_10=subs(goal_para,value_1,-1);
f_11=subs(f_10,findsym(sigema_para_h_temp),H_after);
f_12=subs(f_11,findsym(sigema_para_a_temp),A);
f_13=matlabFunction(f_12,'Vars',{[P1,P2]});
powervec_para_ini=[powervec_cf*0.3,powervec_cf*0.8];
f1_A=(-1)*ones(1,L);
f1_b=(-1)*powervec_cf*L;
f1_ub=powervec_cf*L*ones(1,L);
f1_lb=zeros(1,L);
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
[P,fval,exitflag,output]=fmincon(f_13,powervec_para_ini,f1_A,f1_b,[],[],f1_lb,f1_ub,[],options);

%�õ�P֮��ĵ���
P_para=diag(P);
%���¼���sigema_para
sigema_para=zeros(M,1);
for m=1:M
    sigema_para(m,1)=A(:,m)'*pinv(pinv(P_para)+transpose(H_after(m,:))*conj(H_after(m,:)))*A(:,m);
end
%ֻ���н�������
%��sigema������������
[sigema_para,sigema_para_ind]=sort(real((sigema_para)));
sigema_para_final=sigema_para(L,1);
%���¼���rate
rate_para=zeros(1,L);
for l=1:L
    rate_para(1,l)=max(0,log2(real(P_para(l,l)/sigema_para_final)));
end
per_rate_para(count_h,count_p)=max(rate_para);
sum_rate_para(count_h,count_p)=sum(rate_para);


%��ʼ���㾭�������Ż���succ
%���������
P_succ=P_para;
%��A������sigema�����˳���������
A_trans=transpose(A);
A_para=A_trans(sigema_para_ind,:);
A_para=A_para(1:L,:);
%����H
H_para=H_after(sigema_para_ind,:);
H_para=H_para(1:L,:);

%succ_rate_2,ȫ����
%ȫ���еķ���
%�ж��������з�ʽ
per_temp=perms(1:L);
[num1,~]=size(per_temp);
rate_succ_2=zeros(L,num1);
sigema_succ_2=zeros(L,num1);
for n=1:num1
    A_succ_2=A_para(:,per_temp(n,:));
    H_succ_2=H_para(:,per_temp(n,:));
    powervec_succ_2=P(1,per_temp(n,:));
    P_succ_2=diag(powervec_succ_2);
    %��ʼ����rate
    for m=1:L
        if m==1
            sigema_succ_2(1,n)=A_succ_2(1,:)*pinv(pinv(P_succ_2)+transpose(H_succ_2(1,:))*conj(H_succ_2(1,:)))*A_succ_2(1,:)';
        else
            ins_chol_1=pinv(pinv(P_succ_2)+transpose(H_succ_2(m,:))*conj(H_succ_2(m,:)));
            ins_chol_2=ins_chol_1-ins_chol_1*A_succ_2(1:(m-1),:)'*pinv(A_succ_2(1:(m-1),:)')*pinv(ins_chol_1)*pinv(A_succ_2(1:(m-1),:))*A_succ_2(1:(m-1),:)*ins_chol_1;
            sigema_succ_2(m,n)=abs(real(A_succ_2(m,:)*ins_chol_2*A_succ_2(m,:)'));
        end
    end
    [sigema_succ_2_permu,~]=sort(sigema_succ_2(:,n));
    temp_L=per_temp(n,:);
    for count_temp_L=1:L
        rate_succ_2(temp_L(1,count_temp_L),n)=max(0,real(log2(powervec_succ_2(1,count_temp_L)/sigema_succ_2_permu(count_temp_L,1))));
    end
end
sum_rate_succ_2_temp=sum(rate_succ_2);
[max_value,max_pos]=max(sum_rate_succ_2_temp);
per_rate_succ_2(count_h,count_p)=max(rate_succ_2(:,max_pos));
sum_rate_succ_2(count_h,count_p)=max_value;

%succ_rate_3
%��ÿһ�о��ж��������������ж���
%���н�������(ȷ���û���˳��
rate_succ_3=zeros(1,L);
sigema_succ_3=zeros(L,1);
%���ȶ�BETAA��������
ind_beta_1=zeros(M,L);
ind_beta_user_1=zeros(M,L);
sort_ind_beta_1=zeros(1,L);
for m=1:L
    [~,ind_beta_1(m,:)]=sort(BETAA(m,:),'descend');
    for l=1:L
        [~,ind_beta_user_1(m,l)]=find(ind_beta_1(m,:)==l);
    end
    [~,ind_user_1(m,:)]=sort(ind_beta_user_1(m,:));
    %����˳��ı�A��H��P
    A_succ_3=A_para(:,ind_user_1(m,:));
    H_succ_3=H_para(:,ind_user_1(m,:));
    powervec_succ_3=P(1,ind_user_1(m,:));
    P_succ_3=diag(powervec_succ_3);
    %��ʼ����sigema_succ_3
    if m==1
        sigema_succ_3(m,1)=A_succ_3(1,:)*pinv(pinv(P_succ_3)+transpose(H_succ_3(1,:))*conj(H_succ_3(1,:)))*A_succ_3(1,:)';
    else
        ins_chol_1=pinv(pinv(P_succ_3)+transpose(H_succ_3(m,:))*conj(H_succ_3(m,:)));
        ins_chol_2=ins_chol_1-ins_chol_1*A_succ_3(1:(m-1),:)'*pinv(A_succ_3(1:(m-1),:)')*pinv(ins_chol_1)*pinv(A_succ_3(1:(m-1),:))*A_succ_3(1:(m-1),:)*ins_chol_1;
        sigema_succ_3(m,1)=abs(real(A_succ_3(m,:)*ins_chol_2*A_succ_3(m,:)'));
    end
end%û����
%��ʼ����rate_succ_3
for l=1:L
    rate_succ_3(1,ind_user_1(1,l))=max(0,real(log2(P_succ_3(l,l)/max(sigema_succ_3(1:l,1)))));
end
per_rate_succ_3(count_h,count_p)=max(rate_succ_3);
sum_rate_succ_3(count_h,count_p)=sum(rate_succ_3);

%succ_rate_4
%���н�������(ȷ���û���˳��
rate_succ_4=zeros(1,L);
sigema_succ_4=zeros(L,1);
%���ȶ�BETAA��������
ind_beta_2=zeros(M,L);
ind_beta_user_2=zeros(M,L);
sort_ind_beta_2=zeros(1,L);
for m=1:M
    [~,ind_beta_2(m,:)]=sort(BETAA(m,:));
    for l=1:L
        [~,ind_beta_user_2(m,l)]=find(ind_beta_2(m,:)==l);
    end
end
%�����ÿ���û���������
for l=1:L
    sort_ind_beta_2(1,l)=sum(ind_beta_user_2(:,l));
end
%���û���������
[~,ind_user_2]=sort(sort_ind_beta_2);
%����˳��ı�A��H��P
A_succ_4=A_para(:,ind_user_2);
H_succ_4=H_para(:,ind_user_2);
powervec_succ_4=P(1,ind_user_2);
P_succ_4=diag(powervec_succ_4);
%��ʼ����sigema_succ_1
for m=1:L
    if m==1
        sigema_succ_4(m,1)=A_succ_4(1,:)*pinv(pinv(P_succ_4)+transpose(H_succ_4(1,:))*conj(H_succ_4(1,:)))*A_succ_4(1,:)';
    else
        ins_chol_1=pinv(pinv(P_succ_4)+transpose(H_succ_4(m,:))*conj(H_succ_4(m,:)));
        ins_chol_2=ins_chol_1-ins_chol_1*A_succ_4(1:(m-1),:)'*pinv(A_succ_4(1:(m-1),:)')*pinv(ins_chol_1)*pinv(A_succ_4(1:(m-1),:))*A_succ_4(1:(m-1),:)*ins_chol_1;
        sigema_succ_4(m,1)=abs(real(A_succ_4(m,:)*ins_chol_2*A_succ_4(m,:)'));
    end
end
%��ʼ����rate_succ_1
for l=1:L
    rate_succ_4(1,ind_user_2(1,l))=max(0,real(log2(P_succ_4(l,l)/max(sigema_succ_4(1:l,1)))));
end
per_rate_succ_4(count_h,count_p)=max(rate_succ_4);
sum_rate_succ_4(count_h,count_p)=sum(rate_succ_4);


count_p
end
count_h
end
% for count_p2=1:x_value
%     [num_per_cf,~]=size(find(per_rate_cf(:,count_p2)~=0));
%     [num_per_para,~]=size(find(per_rate_para(:,count_p2)~=0));
%     [num_per_succ_1,~]=size(find(per_rate_succ_1(:,count_p2)~=0));
%     [num_per_succ_2,~]=size(find(per_rate_succ_2(:,count_p2)~=0));
%     [num_per_succ_3,~]=size(find(per_rate_succ_3(:,count_p2)~=0));
%     [num_per_succ_4,~]=size(find(per_rate_succ_3(:,count_p2)~=0));
% 
%     [num_sum_cf,~]=size(find(sum_rate_cf(:,count_p2)~=0));
%     [num_sum_para,~]=size(find(sum_rate_para(:,count_p2)~=0));
%     [num_sum_succ_1,~]=size(find(sum_rate_succ_1(:,count_p2)~=0));
%     [num_sum_succ_2,~]=size(find(sum_rate_succ_2(:,count_p2)~=0));
%     [num_sum_succ_3,~]=size(find(sum_rate_succ_3(:,count_p2)~=0));
%     [num_sum_succ_4,~]=size(find(sum_rate_succ_4(:,count_p2)~=0));
% 
%     final_per_rate_cf(1,count_p2)=sum(per_rate_cf(:,count_p2))/num_per_cf;
%     final_per_rate_para(1,count_p2)=sum(per_rate_para(:,count_p2))/num_per_para;
%     final_per_rate_succ_1(1,count_p2)=sum(per_rate_succ_1(:,count_p2))/num_per_succ_1;
%     final_per_rate_succ_2(1,count_p2)=sum(per_rate_succ_2(:,count_p2))/num_per_succ_2;
%     final_per_rate_succ_3(1,count_p2)=sum(per_rate_succ_3(:,count_p2))/num_per_succ_3;
%     final_per_rate_succ_4(1,count_p2)=sum(per_rate_succ_4(:,count_p2))/num_per_succ_4;
% 
%     final_sum_rate_cf(1,count_p2)=sum(sum_rate_cf(:,count_p2))/num_sum_cf;
%     final_sum_rate_para(1,count_p2)=sum(sum_rate_para(:,count_p2))/num_sum_para;
%     final_sum_rate_succ_1(1,count_p2)=sum(sum_rate_succ_1(:,count_p2))/num_sum_succ_1;
%     final_sum_rate_succ_2(1,count_p2)=sum(sum_rate_succ_2(:,count_p2))/num_sum_succ_2;
%     final_sum_rate_succ_3(1,count_p2)=sum(sum_rate_succ_3(:,count_p2))/num_sum_succ_3;
%     final_sum_rate_succ_4(1,count_p2)=sum(sum_rate_succ_4(:,count_p2))/num_sum_succ_4;
% end
result_per=[mean(per_rate_cf);mean(per_rate_para);mean(per_rate_succ_1);mean(per_rate_succ_2);mean(per_rate_succ_3);mean(per_rate_succ_4)];
result_per=real(result_per);
result_sum=[mean(sum_rate_cf);mean(sum_rate_para);mean(sum_rate_succ_1);mean(sum_rate_succ_2);mean(sum_rate_succ_3);mean(sum_rate_succ_4)];
result_sum=real(result_sum);