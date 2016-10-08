clear all
clc
load Exchange_Data
load BarbCont2016.mat
YC=Data.Adj_Returns;
XC=Data.Exchange_Rates;
DC=Data.Date;

%%
Ntrees=100;
Ncomps=length(eq.permnos);
t_l=250;
h_l=20;
%%
start_index=1;
end_index=length(DC);

SL=end_index-start_index+1;
Y=YC(start_index:end_index,:);
for iy=1:size(Y,2)
    %EY(:,iy)=EWMA(Y(:,iy),'com',120,'min_periods',60);
    for jy=11:size(Y,1)
    EY(jy,iy)=mean(Y(jy-10:jy-1,iy));
    end
end
D=DC(start_index:end_index);
X=XC(start_index:end_index,:);
Yhat_Mat=NaN(size(Y));
time_intervals=t_l+1:h_l:SL-h_l;
Feat_Imp=NaN(Ncomps,length(time_intervals),size(XC,2)+1);
for k=1:length(time_intervals)
    k
    t_ind=time_intervals(k);
    sprintf('X = %d : %d ',t_ind-t_l,t_ind-1)
    sprintf('Y = %d : %d ',t_ind-t_l+1,t_ind)
    sprintf('XP = %d : %d ',t_ind,t_ind+h_l-1)
    sprintf('YP = %d : %d \n',t_ind+1,t_ind+h_l)
    temp_X=X(t_ind-t_l:t_ind-1,:);
    for comp_index=1035%1:Ncomps
        temp_Y=Y(t_ind-t_l+1:t_ind,comp_index);
        temp_EY=EY(t_ind-t_l:t_ind-1,comp_index);
        INDICES=~isnan(temp_Y);
        if(sum(INDICES)>t_l*0.2)
            input_Y=temp_Y(INDICES);
            input_X=temp_X(INDICES,:);
            jj=1;
            for ii=1:length(INDICES)
                if (INDICES(ii)==1)
                    input_X(jj,size(XC,2)+1)=temp_EY(ii);
                    jj=jj+1;
                end
            end
            T=TreeBagger(Ntrees,input_X,input_Y,'method','regression','oobvarimp','on');
            Feat_Imp(comp_index,k,:)=T.OOBPermutedVarDeltaError;
            XP=X(t_ind:t_ind+h_l-1,:);
            for ii=1:h_l
                XP(ii,size(XC,2)+1)=EY(ii+t_ind-1,comp_index);
            end
            [Yhat,stdevs] =predict(T,XP);
            Yhat_Mat(t_ind+1:t_ind+h_l,comp_index)=Yhat(1:h_l);            
        end
    end
end
Results.R=Y;
Results.Rhat=Yhat_Mat;
Results.Date=D;
Results.Feat_Imp=Feat_Imp;
save('Results','Results')
save bellairs2
Performance=EvalPerf(Y(:,comp_index),Yhat_Mat(:,comp_index),D);