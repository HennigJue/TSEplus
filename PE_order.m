function [PEorder,PEindex,NS,IPS] = PE_order(PEtype,NY,Necho,k0,PE_steps,HFfac,Nref,GRfac,PEcont)
if(nargin<9), PEcont=0; end

if(~isempty(PE_steps))
    PE_steps=sort(PE_steps,'ascend');
    PEtype='linear';
    disp('PE-ype changed to linear encoding');
end

if(isempty(PE_steps))
    PE_first=floor(1+HFfac*floor(NY/2));
    PE_ref=floor(NY/2-Nref/2);
    GRPElow=[PE_ref:-GRfac:PE_first];
    PE_steps=[GRPElow(end:-1:1) PE_ref+1:PE_ref+Nref-1 PE_ref+Nref:GRfac:NY];
    nPE=length(PE_steps);
    NS=floor(nPE/Necho);
    if(NS==0), NS=1; Necho=length(PE_steps);
    else
        PE_steps=PE_steps(1:NS*Necho);
    end
else
    PE_steps=PE_steps+floor(NY/2);
end
nPE=length(PE_steps);
NS=floor(nPE/Necho);
if(NS==0), NS=1; Necho=length(PE_steps);
else
    PE_steps=PE_steps(1:NS*Necho);
end
ind=find(PE_steps<NY/2);
nPE0=ind(end)+1;
PEindex=PE_steps;
if(NS==1)
    if(strcmp(PEtype,'faise')),
        disp('PE-mode changed to centric encoding');
        PEtype='centric';
    end
end


switch PEtype
    case 'linear'
        %PE0=ind(end)+1;
        PEorder=reshape(PE_steps,[NS Necho])'-floor(NY/2);
        [ke,ks]=find(floor(PEorder)==0);
        PEorder=circshift(PEorder,[k0-ke(1) 0]);

    case {'centric','centric_s'}
        PEcont=0;
        kn=PE_steps-floor(NY/2);
        ind=find(kn<0);
        kn(ind)=-kn(ind)+0.5;
        kn=sort(kn);
        temp=kn-floor(kn);
        ind1=find(temp>0);
        kn(ind1)=-kn(ind1)+0.5;
        kn=kn-min(kn)+1;
        kn=kn-nPE0;
        knorder=reshape(kn,[NS Necho]);
        if(strcmp(PEtype,'centric_s'))
            PEorder=circshift(knorder,[0 k0]);
        else
            PEorder=knorder-knorder(floor(NS/2+1),k0);
        end
        PEorder=PEorder';
        ind=find(PEorder<-NY/2);




    case 'faise'
        if(floor(NS/2)==NS/2)
            kn=0*PE_steps;
            ind=0;
            kn(1)=ind(1);
            for k=1:floor(length(PE_steps)/2)
                kn(2*(k-1)+2)=ind+k;
                kn(2*(k-1)+3)=ind-k;
            end
            NS=floor(length(PE_steps)/Necho);
            kn=kn(1:NS*Necho);
            kn=kn-min(kn)+1;
            knorder=reshape(kn,[NS Necho]);
            [ke,ks]=find(knorder'==nPE0);
            knorder=knorder+NS/2*(k0-ke);
            knorder=mod(knorder-1,length(PE_steps))+1;
            PEorder=PE_steps(knorder);
            PEorder=PEorder'-floor(NY/2);
            [ke,ks]=find(knorder'==nPE0);

        else
            NYF=nPE;
            IPS=round(NS/4)+round((k0-1)*NS/2);
            for m=1:floor(NS/2)
                for n=1:(NYF/NS)
                    A(m,n)=mod(NYF/2-n*NS/2+m-1+IPS,NYF)+1;
                end
            end
            for m=floor(NS/2)+1:NS
                for n=1:(NYF/NS)
                    A(m,n)=mod(NYF/2+(n-2)*NS/2+m-1+IPS,NYF)+1;
                end
            end
            for m=1:NS
                for n=1:floor(NYF/NS)
                    k(m,n)=(A(m,n)-(NYF+1)/2);
                end
            end
            ks=k(:);
            [kss,indk]=sort(ks-min(ks)+1);
            PE_new(indk)=PE_steps;
            PE_new1=PE_new-PE_new(1)+ks(1);
            PEorder=reshape(PE_new1,[NS Necho]);
            PEorder=PEorder';
        end

end
if(PEcont==1)
    kk=0;
    dPEorder=diff(PEorder(k0:end,:));
    [inde,inds]=find(abs(dPEorder)>NY/2);
    for k=1:length(inds)
        PEorder(inde(k)+k0:end,inds(k))=PEorder(inde(k)+k0:end,inds(k))+sign(dPEorder(inde(k),inds(k)))*dPEorder(inde(k),inds(k))+dPEorder(inde(k)-1,inds(k));
    end
else
    %PEorder=mod(PEorder,NY)-round(NY/2);
end
% if(PEcont==0)
% 
% else
%     PEorder(ind)=PEorder(ind)+NY;
% end