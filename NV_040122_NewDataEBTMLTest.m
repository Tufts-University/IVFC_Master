%% New Data Testing
% Use this code to evaluate if how well new data performs on the trained
% EBT model.
A=readcell('reference.csv');
%% Remove Training and Test(Val) Data
A(18:71,:)=[];
A(1,:)=[];
%% Load trained model
fname=A(:,1);
ftype=A(:,2);
for i=1:38
    fname{i}=fname{i}(2:end-1);
    ftype{i}=ftype{i}(2:end-1);
end
% load('GBEnsemble_model89.mat')
% filenum=filenum';
filenum=nonzeros(filenum(:));
%% Formatting Loop
label_store_TP=[];
r_store_TP=[];
peak_value_store_TP=[];
label_store_FP=[];
r_store_FP=[];
peak_value_store_FP=[];

for lm=[1:38]
    k=filenum(lm);
    disp(['Currently on Day # ',num2str(round(lm/3)),' of ', num2str(size...
        (fname,1)./3)])
    [range_data,label,peak_values]=MLformating6(fname2{lm},ftype{lm});
    lbead=find(label==3); %Remove any beads from data
    label(lbead)=[];
    range_data(lbead,:)=[];
    peak_values(lbead,:)=[];

    lTP=find(label==1); %Remove TPs from data
    label(lTP)=[];
    range_data(lTP,:)=[];
    peak_values(lTP,:)=[];

    peak_values(:,22)=lm;
    X=[range_data,label];
    label_store_FP=[label_store_FP;X(:,109)];
    r_store_FP=[r_store_FP;X(:,1:81)];
    peak_value_store_FP=[peak_value_store_FP;peak_values];

    [range_data,label,peak_values]=MLformating6(fname{lm},ftype{lm});
    lbead=find(label==3);
    label(lbead)=[];
    range_data(lbead,:)=[];
    peak_values(lbead,:)=[];

    peak_values(:,22)=lm;
    X=[range_data,label];
    label_store_TP=[label_store_TP;X(:,109)];
    r_store_TP=[r_store_TP;X(:,1:81)];
    peak_value_store_TP=[peak_value_store_TP;peak_values];
end
X_TP_store=[r_store_TP,label_store_TP];
X_FP_store=[r_store_FP,label_store_FP];
%%
PV_testing=[peak_value_store_TP;peak_value_store_FP];
X_testing=[X_TP_store;X_FP_store];
%% Testing Loop
X=X_testing;
PV=PV_testing;
tAcc=[];
pur_store=[];
b_store = cell(50, 1);
Peak_values_deleted=[];
for iii=1:50
    disp(['Testing: Running Model # ',num2str(iii),' of 50'])
    X_test=X(:,1:81);
    label=X(:,82);
    yfit = GBEnsemble{1}(iii).predictFcn(X_test);
    a=label-yfit;
    b=[length(find(label(a==0)==1)),length(find(a==1));length(find(a==-1)),length(find(label(a==0)==2))];
    Sens(iii)=(b(1,1)./(b(1,1)+b(2,1)));
    Spec(iii)=(b(2,2)./(b(2,2)+b(1,2)));
    acc=(b(1,1)+b(2,2))./(sum(sum(b)));
    b_store{iii}=b;
    tAcc=[tAcc;acc];
    pur_store=[pur_store;(b(1,1)./(b(1,1)+b(1,2)))];
    X=[X_test(find(yfit==1),:),label(find(yfit==1))];
    Peak_values_deleted=[Peak_values_deleted;PV(find(yfit==2),:)];
    Peak_values_saved=PV(find(yfit==1),:);
    PV=PV(find(yfit==1),:);
end
%%
mat_pur=[];
mat_sens=[];
mat_spec=[];
mat_acc=[];
for lm=1:50
    net_pos=0;
    net_neg=0;
    for k=1:lm
        b22=b_store{k}(2,2);
        b12=b_store{k}(1,2);
        b11=b_store{k}(1,1);
        b21=b_store{k}(2,1);
        net_pos=b21+net_pos;
        net_neg=b22+net_neg;
    end
    sen=b11/(b11+net_pos);
    mat_sens=[mat_sens;sen];
    spe=net_neg/(b12+net_neg);
    mat_spec=[mat_spec;spe];
    pur=b11/(b11+b12);
    mat_pur=[mat_pur;pur];
    acc=(b(1,1)+net_neg)./(net_neg+net_pos+b(1,1)+b(1,2));
    mat_acc=[mat_acc;acc];
end
total_mat=[mat_pur,mat_spec,mat_sens,mat_acc];
%% ML Formatting Code
function [range_data,label,peak_values]=MLformating6(fname,ftype)
Wn=[50 6E3]./(60E3/2);%Cutoff frequencies divided by Nyquist frequency
[b,a]=butter(2,Wn);
load(fname)
fileN=ftype;
Spec5=[7.10E+01	7.07E+01	7.03E+01];
peak_values(find(peak_values(:,9)<20),:)=[];
peak_values=sortrows(peak_values,[21,6,7]);
chunk=unique(peak_values(:,6));
range_data=zeros(size(peak_values,1),108);
label=2.*ones(size(peak_values,1),1);
count=1;
for i=1:max(chunk)
    l1=find(peak_values(:,6)==i);
    load([fileN,'_',num2str(i),'_raw.mat'])
    load([fileN,'_sigmas']);


    %     M(:,1)=M(:,1)./abs(Spec5(1));
    %     M(:,2)=M(:,2)./abs(Spec5(2));
    %     M(:,3)=M(:,3)./abs(Spec5(3));
    %     M(:,4)=M(:,4);
    %     M(:,5)=M(:,5);
    avg_col=mean(M);
    std_col=std(M);
    %     M(:,1)=(M(:,1))./mean(M(M(:,1)>avg_col(:,1)+5.*std_col(:,1),1));
    %     M(:,2)=(M(:,2))./mean(M(M(:,2)>avg_col(:,2)+5.*std_col(:,2),2));
    %     M(:,3)=(M(:,3))./mean(M(M(:,3)>avg_col(:,3)+5.*std_col(:,3),3));
    M(:,1)=(M(:,1)-median(M(:,1)))./std(M(:,1));
    M(:,2)=(M(:,2)-median(M(:,2)))./std(M(:,2));
    M(:,3)=(M(:,3)-median(M(:,3)))./std(M(:,3));
    M_filt(:,4)=filtfilt(b,a,M(:,4));
    M(:,4)=(M_filt(:,4)-median(M_filt(:,4)))./std(M_filt(:,4));
    M_tot=M(:,1)+M(:,2)+M(:,3);
    %
    %     M(:,1)=(M(:,1)-mean(M(:,1)));
    %     M(:,2)=(M(:,2)-mean(M(:,2)));
    %     M(:,3)=(M(:,3)-mean(M(:,3)));
    %     M(:,4)=(M(:,4)-mean(M(:,4)));
    %     M(:,5)=(M(:,5)-mean(M(:,5)));
    %     % Normalizing by standard deviation
    %
    %     SN_405=(M(:,1)-mean(M(:,1)))./sigmas_final(1);
    %     SN_488=(M(:,2)-mean(M(:,2)))./sigmas_final(2);
    %     SN_633=(M(:,3)-mean(M(:,3)))./sigmas_final(3);
    %     SN_Red=(M(:,4)-mean(M(:,4)));
    %     SN_Green=(M(:,5)-mean(M(:,5)));
    %     clear M
    %     M=[SN_405,SN_488,SN_633,SN_Red,SN_Green];

    disp(['Chunk # ',num2str(i),' of ',num2str(max(chunk))])
    for j=1:size(l1,1)
        loc=peak_values(l1(j),7);
        M_range=[M(loc-13:loc+13,1)',M(loc-13:loc+13,2)',...
            M(loc-13:loc+13,3)',M(loc-13:loc+13,4)'];
        if peak_values(l1(j),4)<0.3&peak_values(l1(j),5)>0.15
            %if peak_values(l1(j),5)>0.2
            label(count)=1;
        elseif peak_values(l1(j),4)>0.3
            label(count)=3;
        end
        range_data(count,:)=M_range;
        count=count+1;
    end
    clear M_filt
end
end