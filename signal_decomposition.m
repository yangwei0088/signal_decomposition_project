function results=signal_decomposition(ml,mr,k)

%   Parameters:
%       -----------------
%       ml: parameter for penalizing within period smoothness
%       mr: parameter for penalizing between period smoothness
%       k: parameter for adjusting the low-rank structure of the mean
%       signal
%       -----------------
%  Output:
%       -----------------
%       results   :  a struct containing the reconstructed mean,
%       anomaly, and original signal

%default parameter values;
ml=100000;
mr=6000000;
k=3;

all_files_numbers=26:26+19;
current=26+19;
for i=1:8
    all_files_numbers=[all_files_numbers current+6:current+6+19];
    current=current+6+19;
end

addpath('functions');
folder_path=strcat('/input/','input1');
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


ratios=cell([],1);
for iter=all_files_numbers
    data = readtable(filenames{iter},'Delimiter', ',');
    power=data.Power;
    poa=data.POA;
    fault=data.anomaly;
    ratio=power./poa;
    ratio(ratio<0)=0;
    ratio(isnan(ratio))=0;
    ratios{iter}=ratio;
    M0=reshape(ratio(1:10080),1440,7);
    M0(isnan(M0))=0;
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

    % minimize U and V 
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
                 Z=A+omega*V';
                 omega=max(omega+alpha*(-omega*(VTV)-A*V),0);
                lagrangian_objective1(i)=0.5*sum((omega*V.').^2,'all')-sum(omega.*((A+omega*V.')*V),'all');
            end
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
            ro=0.000001;   
            V  = sylvester(U.'*U+ro*eye(size(U,2)),mr*(D2.'*D2),U.'*B+ro*Z-ro*hat);
            A=V+hat;
            lambda=rand(size(M0,1),size(M0,2));
            lagrangian_objective=zeros(50,1);
            for i=1:50
                %alpha=alpha*0.2;
                Z=A+ U.'*lambda;
                lambda=max(lambda+alpha1*(-(UUT)*lambda-U*A),0);
                %lagrangian_objective(i)=0.5*sum((U.'*lambda).^2,'all')-sum(lambda.*(U*(A+U.'*lambda)),'all');
            end
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
    results.anomaly_detected{iter}=temp1(:);
    temp=U*V;
    temp=max(0,temp);
    results.reconstruct{iter}=temp(:);
    results.original{iter}=M0(:); 
    results.anomaly_true{iter}=fault;
    
end
end

