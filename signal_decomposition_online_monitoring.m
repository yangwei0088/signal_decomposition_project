
% default values
ml=100000;
mr=6000000;
k=3;

addpath('functions');
folder_path=strcat('input/','input2');
disp(folder_path);
a=dir(folder_path);
foldernames=cell([],1);
count=1;
for i=1:length(a)
    if a(i).name(1)~='.'
        foldernames{count}=a(i).name;
        count=count+1;
    end
end
foldernames=natsort(foldernames);
count=1;
filenames=cell([],1);
for f=1:length(foldernames)
    newpath=strcat(folder_path,'/',foldernames(f));
    a=dir(newpath{1});
    for i=1:length(a)
        if a(i).name(1)~='.'
            temp=strcat(folder_path,'/',foldernames(f),'/',a(i).name);
            filenames{count}=temp{1};
            count=count+1;
        end
    end
end
%%
load('normal6.mat');
%%
threshold=31.2;
monitor_count=0;
total_count2=0;
total_count3=0;
count2=0;
count3=0;

for iter=1:length(filenames)
    disp([iter,length(filenames)]);
    data = readtable(filenames{iter},'Delimiter', ',');
    power=data.Power;
    poa=data.POA;
    fault=data.anomaly;
    ratio=power./poa;
    ratio(ratio<0)=0;
    ratio(isnan(ratio))=0;
    for day_iter=1:7
        power1=power(1+(day_iter-1)*1440:1440*day_iter);
        poa1=poa(1+(day_iter-1)*1440:1440*day_iter);
        %poa1(isnan(poa1))=0;
        start_idx=find(power1,1,'first');
        end_idx=find(power1,1,'last');
        temp1=matches(fault(1440*(day_iter-1)+1:1440*(day_iter)),'True');

        if (sum(temp1)>=1) || (sum(poa1(start_idx+100:end_idx-100)<3)==0)
            M0=[normal6;ratio(1440*(day_iter-1)+1:1440*day_iter)];
            M0=reshape(M0,1440,7);
            M0(isnan(M0))=0;
            monitor_count=monitor_count+1;
            disp(monitor_count);

            [U0,S,V]=svd(M0);
            U=U0(:,1:k)*S(1:k,1:k);
            V=S(1:k,1:k)*V(1:k,:);
            Im=eye(k);
            In=eye(size(M0,1));
            Imn=eye(k*(size(M0,1)));
            D=zeros(size(M0,1)-2,size(M0,1));
            hat=zeros(size(M0,1),k);
            omega=zeros(size(M0,1),size(M0,2));
            omega_plus=omega+1;
            for i=1:(size(M0,1)-2)
                for j=i:i+2
        
                    if j==i || j==i+2
                        D(i,j)=1;
                    else
                        D(i,j)=-2;
                    end 
                end
            end
        
            % build D2
             D2=zeros(size(M0,2)-2,size(M0,2));
            for i=1:(size(M0,2)-2)
                for j=i:i+2
                    if j==i || j==i+2
                        D2(i,j)=1;
                    else
                        D2(i,j)=-2;
                    end 
                end
            end
        
             overall_objective=zeros([],1);
             n=size(M0,1);
             order=1;
             knots=[ones(1,order) linspace(1,n,40) n*ones(1,order)];
             nknots=length(knots)-order;
             kspline=spmak(knots,eye(nknots));
             H1=spval(kspline,1:n)';
        
             n=size(M0,2);
             knots=[ones(1,order) linspace(1,n,10) n*ones(1,order)];
             nknots=length(knots)-order;
             kspline=spmak(knots,eye(nknots));
             H2=spval(kspline,1:n)';
             theta=zeros(size(H1,2),size(H2,2));
        
             %
             alpha=1*10^(-9);
             alpha1=0.1;
            for q=1:5
                B=M0-H1*theta*H2.';
                disp(q); 
                VTV=V'*V; 
                hat=rand(size(M0,1),k);
                Z=rand(size(M0,1),k);
                result = sum(isnan(V(:)));
                if result==0
                    L=max(eig(VTV));
                    alpha=1/L;
                end
                %alpha=;
                for t=1:10
                    
                    ro=0.000001;
                    omega=rand(size(M0,1),size(M0,2));
                    U = sylvester(ml*(D.'*D)+ro*eye(size(D,2)),V*V.',B*V.'+ro*Z-ro*hat);
                    A=U+hat;
                    lagrangian_objective1=zeros(20,1);
                    for i=1:50
                         
                         omega=max(omega+alpha*(-omega*(VTV)-A*V),0);
                        lagrangian_objective1(i)=0.5*sum((omega*V.').^2,'all')-sum(omega.*((A+omega*V.')*V),'all');
                    end
                    Z=A+omega*V';
                    hat=hat+U-Z;
                end
                Z=zeros(k,size(M0,2));
                hat=rand(k,size(M0,2));
                UUT=U*U';
                %L=max(eig(UUT));
                result = sum(isnan(U(:)));
                if result==0
                    L=max(eig(UUT));
                    alpha1=1/L;
                end
                for t=1:10
                    %alpha1=alpha1+100;
                    ro=0.000001;
                   
                    V  = sylvester(U.'*U+ro*eye(size(U,2)),mr*(D2.'*D2),U.'*B+ro*Z-ro*hat);
                    A=V+hat;
                    lambda=rand(size(M0,1),size(M0,2));
                    lagrangian_objective=zeros(50,1);
                    for i=1:50
                        %alpha=alpha*0.2;
                        
                        lambda=max(lambda+alpha1*(-(UUT)*lambda-U*A),0);
                        %lagrangian_objective(i)=0.5*sum((U.'*lambda).^2,'all')-sum(lambda.*(U*(A+U.'*lambda)),'all');
                    end
                    Z=A+ U.'*lambda;
                   %temp1=U*Z;
                   %disp(sum(temp1(temp1<0)));
                    %Z=max(A,0);
                    hat=hat+V-Z;
                end
                %disp(0.5*sum((U*V-B).^2,'all')+0.5*ml*sum((D*U).^2,'all')+0.5*mr*sum((D2*V.').^2,'all'));
                %overall_objective(n)=0.5*sum((U*V-B).^2,'all')+0.5*ml*sum((D*U).^2,'all')+0.5*mr*sum((D2*V.').^2,'all');
        
             Y=M0-U*V;
             %
            penalty_theta=1*10^(1.2);
            alpha1=3*10^(-2);
            alpha2=3*10^(-2);
            overall_objective=zeros([],1);
            ro=10^(-2);
             % constrained least square with lasso penalty step 1 
             Z=zeros(size(H1,2),size(H2,2));
             u=zeros(size(H1,2),size(H2,2));
             thr=penalty_theta*alpha1/2;
             for i=1:20
                 theta_old=theta;
                 theta=zeros(size(H1,2),size(H2,2));
                 objective1=zeros([],1);
                 for t=1:500
                    %alpha1=alpha1*0.99;
                    thr=penalty_theta*alpha1/2;
                    temp=H1.'*Y*H2-(H1.'*H1)*theta*(H2.'*H2)-ro*(theta-Z+u);
                    temp1=theta+alpha1*temp;
                    theta=wthresh(temp1,'s',thr);
                    level = graythresh(temp1);
                    penalty_theta=level/alpha1;
                    % theta=wthresh(theta,'s',thr);
                    objective1(t)=sum((Y-H1*theta*H2.').^2,'all')+penalty_theta*norm(theta,1)+ro/2*sum((theta-Z+u).^2,'all');
                 end
                 theta_new=theta;
                 A=theta+u;
                 objective2=zeros([],1);
                 lambda=zeros(size(H1,1),size(H2,1));
        
            %     cvx_end
                 for j=1:100
                     lambda=max(lambda+alpha2*(-(H1*H1.')*lambda*(H2*H2.')+H1*A*H2.'),0);
                     objective2(j)=0.5*sum((H1.'*lambda*H2).^2,'all')+sum(lambda.*(H1*A*H2.'),'all')-sum(lambda.*((H1*H1.')*lambda*(H2*H2.')),'all');
                 end
                 Z=A-H1.'*lambda*H2;
                 u=u+theta-Z;
                 %overall_objective(i)=0.5*sum((Y-H1*theta*H2.').^2,'all')+penalty_theta*norm(theta,1);
                 %disp(overall_objective(i));
                 %disp(sum(theta_new-theta_old,'all'));
             end
            end

            temp1=H1*theta*H2.';
            temp1=min(temp1,0);
            results.anomaly_detected{monitor_count}=temp1(:);
            temp=U*V;
            temp=max(0,temp);
            results.reconstruct{monitor_count}=temp(:);
            results.original{monitor_count}=M0(:); 
            results.anomaly_true{monitor_count}=fault;
            
            anomaly_detected=results.anomaly_detected{monitor_count}(1440*6+1:1440*7);
            start_idx=find(anomaly_detected,1,'first');
            end_idx=find(anomaly_detected,1,'last');
            anomaly_detected(1:start_idx)=0;
            anomaly_detected(end_idx:end)=0;
            temp=anomaly_detected;
            temp1=matches(fault(1440*(day_iter-1)+1:1440*(day_iter)),'True');
            location=0;
            temp33=temp;
            poa11=poa1;
            poa11(isnan(poa11))=0;
            start_idx=find(poa11,1,'first');
            end_idx=find(poa11,1,'last');
            temp33(1:start_idx)=0;
            temp33(end_idx:end)=0;
            t=start_idx;
            while t<end_idx
                if temp33(t)<0
                    h=1;
                    while t+h<end_idx && temp33(t+h)<0
                        temp33(t+h)=temp33(t+h)+temp33(t+h-1);
                        h=h+1;
                    end
                    t=t+h-1;
                end
                t=t+1;
            end
            
            if sum(temp1)>=1
                total_count2=total_count2+1;
                if sum(temp33<-threshold)>1
                    count2=count2+1;
                else
                    normal6=[normal6(1+1440:end);M0(:,7)];
                end
                disp('detection');
                disp([count2,total_count2]);
            else 
                start_idx=find(power1,1,'first');
                end_idx=find(power1,1,'last');
                if sum(poa1(start_idx+100:end_idx-100)<3)==0
                    total_count3=total_count3+1;
                    if sum(temp33<-threshold)>1
                        count3=count3+1;
                    else
                        normal6=[normal6(1+1440:end);M0(:,7)];
                    end  
                end
                disp('false alarm');
                disp([count3,total_count3]);
            end
        end
    end
end






