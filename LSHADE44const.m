function [vystup] = LSHADE44const(Dim,Max_Gen,xmin,xmax,fnum,eps_viol)%
rand('state',sum(100*clock));
%binomialni a exponencialni krizeni soutezi ale jen v jedne generaci
global initial_flag

N_init=18*Dim;
N_min=4;
ps=N_init;
p=0.11;
max_velikost_archivu=round(ps*2.6);

 
pesek1=0;pesek2=0;
 

 maxiter=Max_Gen*Dim;
 evals=ps;
 D=Dim;
 H=6;

 h=4;
 n0=2;
 delta=1/(5*h);
 ni=zeros(1,h)+n0;

 Xmin=repmat(xmin,1,D);
 Xmax=repmat(xmax,1,D);
 Xmin=repmat(Xmin,ps,1);
 Xmax=repmat(Xmax,ps,1);
 
 pos=Xmin+(Xmax-Xmin).*rand(ps,D);

 P=zeros(ps,D+11);
 P(:,1:D)=pos;
 initial_flag = 0;
 [P(:,D+1),gf,hf]=CEC2017(pos,fnum);
 [P(:,D+2),P(:,D+3:D+8),P(:,D+9:D+11)]=violation_velke(fnum,hf,gf,eps_viol,ps);  
 
%  [violmin, indexviolmin]=min(P(:,D+2)); %minimalizuju violation
%  pomocna=P(indexviolmin,:);
%  fmin=min(pomocna(:,D+1));

 
 MFpbin=.5*ones(1,H);
 MFpexp=.5*ones(1,H);
 MCRpbin=.5*ones(1,H);%neni zde CR ale pm
 MCRpexp=.5*ones(1,H);
 MFprlbin=.5*ones(1,H);
 MFprlexp=.5*ones(1,H);
 MCRprlbin=.5*ones(1,H);%neni zde CR ale pm
 MCRprlexp=.5*ones(1,H);
 %pro ulozeni prumernych hodnot v pameti (vsechny MCR a vsechny MF)
 %pro ulozeni Fmin-Fmax v 10 etapach
  
 kpbin=1;
 kpexp=1;
 kprlbin=1;
 kprlexp=1;
 
 velarchivu=0;
 A=[];

 while (evals<maxiter) 
    Fpoleexp=-1*ones(1,ps);
    CRpole=-1*ones(1,ps);
    CRpoleexp=-1*ones(1,ps);
    Fpolerl=-1*ones(1,ps);
    Fpolerlexp=-1*ones(1,ps);
    CRpolerl=-1*ones(1,ps);
    CRpolerlexp=-1*ones(1,ps);
    
    strategie=zeros(1,ps);
    
    SCRpbin=[];SFpbin=[];
    SCRpexp=[];SFpexp=[];
    SCRprlbin=[];SFprlbin=[];
    SCRprlexp=[];SFprlexp=[];
    uspesnychpbin=0;
    uspesnychpexp=0;
    uspesnychprlbin=0;
    uspesnychprlexp=0;
    
    poskon=zeros(ps,D);

    deltafcepbin=-1*ones(1,ps);
    deltafcepexp=-1*ones(1,ps);
    deltafceprlbin=-1*ones(1,ps);
    deltafceprlexp=-1*ones(1,ps);
    
    for i=1:ps  %VYTVORENI DALSI GENERACE
        [hh,p_min]=roulete(ni);
        if p_min<delta
            ni=zeros(1,h)+n0;
        end  %reset
        r=nahvyb(H,1) ;%        r=1+fix(H*rand(1));
        switch hh
          case 1 %(CURRENTTORAND/BIN)
           strategie(1,i)=1;
           if MCRpbin(1,r)==-1 
               CR=0;
            else
               CR=MCRpbin(1,r)+ sqrt(0.1)*randn(1);
            end
            if CR>1
                CR=1;
            else if CR<0
                CR=0;
                end
            end
            F=-1;
            while F<=0  
                F=rand*pi-pi/2; 
                F=0.1 * tan(F) + MFpbin(1,r);   
            end
            if F>1
                F=1;
            end
%             p = pmin+ (0.2-pmin) * rand;
            ppoc=round(p*ps);
            if ppoc==0 
                ppoc=1;
            end
            Fpole(1,i)=F;
            CRpole(1,i)=CR;
            y=currenttopbestbin_izrc_cons(A,velarchivu,P(:,1:D),P(:,D+1),P(:,D+2),ppoc,F,CR,i,xmin,xmax);
            poskon(i,:)=y;


          case 2  %(CURRENTTORAND/EXP)
            strategie(1,i)=2;
           if MCRpexp(1,r)==-1 
               CR=0;
            else
               CR=MCRpexp(1,r)+ sqrt(0.1)*randn(1);
            end
           if CR>1
                CR=1;
            else if CR<0
                CR=0;
                end
           end
            F=-1;
            while F<=0  
                F=rand*pi-pi/2; 
                F=0.1 * tan(F) + MFpexp(1,r);   
            end
            if F>1
                F=1;
            end
%             p = pmin+ (0.2-pmin) * rand;
            ppoc=round(p*ps);
            if ppoc==0 
                ppoc=1;
            end
            Fpoleexp(1,i)=F;
            CRpoleexp(1,i)=CR;
            y=currenttopbestexp_izrc_cons(A,velarchivu,P(:,1:D),P(:,D+1),P(:,D+2),ppoc,F,CR,i,xmin,xmax);
            poskon(i,:)=y;

            case 3  %(RANDRL/BIN)
           strategie(1,i)=3;
           if MCRprlbin(1,r)==-1 
               CR=0;
            else
               CR=MCRprlbin(1,r)+ sqrt(0.1)*randn(1);
            end
            if CR>1
                CR=1;
            else if CR<0
                CR=0;
                end
            end
            F=-1;
            while F<=0  
                F=rand*pi-pi/2; 
                F=0.1 * tan(F) + MFprlbin(1,r);   
            end
            if F>1
                F=1;
            end
            
            
            Fpolerl(1,i)=F;
            CRpolerl(1,i)=CR;
            y=derand_RLe_cons(P(:,1:D),P(:,D+1),P(:,D+2),F,CR,i);
            poskon(i,:)=zrcad(y,xmin,xmax);
            
            case 4  %(RANDRL/EXP)
            strategie(1,i)=4;
           if MCRprlexp(1,r)==-1 
               CR=0;
            else
               CR=MCRprlexp(1,r)+ sqrt(0.1)*randn(1);
            end
           if CR>1
                CR=1;
            else if CR<0
                CR=0;
                end
           end
            F=-1;
            while F<=0  
                F=rand*pi-pi/2; 
                F=0.1 * tan(F) + MFprlexp(1,r);   
            end
            if F>1
                F=1;
            end
            Fpolerlexp(1,i)=F;
            CRpolerlexp(1,i)=CR;
            y=derandexp_RLe_cons(P(:,1:D),P(:,D+1),P(:,D+2),F,CR,i);
            poskon(i,:)=zrcad(y,xmin,xmax);
        
        end            
        
    end
   Q=zeros(ps,D+11);
   Q(:,1:D)=poskon;
   initial_flag = 0;
   [Q(:,D+1),gf,hf]=CEC2017(poskon,fnum);
   [Q(:,D+2),Q(:,D+3:D+8),Q(:,D+9:D+11)]=violation_velke(fnum,hf,gf,eps_viol,ps);  
    
% zjisteni, jak jsou na tom prvky Q

   jak=zeros(ps); 
   for i=1:ps
       if  Q(i,D+2)==0 && P(i,D+2)==0
           if Q(i,D+1)<= P(i,D+1)
              % nahrad - y je uspesny
              jak(i)=1; 
           end
       elseif Q(i,D+2) <= P(i,D+2)
              % nahrad
              jak(i)=1; 
       end
   end
    for i=1:ps
       if jak(i)==1 
            switch  strategie(1,i) 
                case 1
                    if Q(i,D+2)==0 && P(i,D+2)==0
                        deltafcepbin(1,i)=P(i,D+1)-Q(i,D+1);
                    else
                        deltafcepbin(1,i)=P(i,D+2)-Q(i,D+2);
                    end
                    uspesnychpbin=uspesnychpbin+1;
                    SCRpbin=[SCRpbin,CRpole(1,i)];
                    SFpbin=[SFpbin,Fpole(1,i)];
                case 2
                    if Q(i,D+2)==0 && P(i,D+2)==0
                        deltafcepexp(1,i)=P(i,D+1)-Q(i,D+1);
                    else
                        deltafcepexp(1,i)=P(i,D+2)-Q(i,D+2);
                    end
                    uspesnychpexp=uspesnychpexp+1;
                    SCRpexp=[SCRpexp,CRpoleexp(1,i)];
                    SFpexp=[SFpexp,Fpoleexp(1,i)];
                case 3 
                    if Q(i,D+2)==0 && P(i,D+2)==0
                        deltafceprlbin(1,i)=P(i,D+1)-Q(i,D+1);
                    else
                        deltafceprlbin(1,i)=P(i,D+2)-Q(i,D+2);
                    end
                    uspesnychprlbin=uspesnychprlbin+1;
                    SCRprlbin=[SCRprlbin,CRpolerl(1,i)];
                    SFprlbin=[SFprlbin,Fpolerl(1,i)];
                    
                case 4
                    if Q(i,D+2)==0 && P(i,D+2)==0
                        deltafceprlexp(1,i)=P(i,D+1)-Q(i,D+1);
                    else
                        deltafceprlexp(1,i)=P(i,D+2)-Q(i,D+2);
                    end
                    uspesnychprlexp=uspesnychprlexp+1;
                    SCRprlexp=[SCRprlexp,CRpolerlexp(1,i)];
                    SFprlexp=[SFprlexp,Fpolerlexp(1,i)];
            end
            if velarchivu < max_velikost_archivu
                A=[A;P(i,1:D)];
                velarchivu=velarchivu+1;
            else
                ktere=nahvyb(velarchivu,1);
                A(ktere,:)=P(i,1:D);
            end
            P(i,:)=Q(i,:);
            ni(strategie(1,i))=ni(strategie(1,i))+1;         % zmena prsti qi
       end
    end
    
    if uspesnychpbin>0
        platne=find(deltafcepbin~=-1);
        delty=deltafcepbin(1,platne);
        suma=sum(delty);
        vahyw=1/suma*delty;
        mSCRpbin=max(SCRpbin);
        if MCRpbin(1,kpbin)==-1  ||  mSCRpbin==0
            MCRpbin(1,kpbin)=-1;
        else    
            MCRpbin(1,kpbin)=sum(vahyw.*SCRpbin);
        end
        meanSFpomjm=vahyw.*SFpbin;
        meanSFpomci=meanSFpomjm.*SFpbin;
        MFpbin(1,kpbin)=sum(meanSFpomci)/sum(meanSFpomjm);
        kpbin=kpbin+1;
        if kpbin>H
            kpbin=1;
        end
    end
    if uspesnychpexp>0
        platne=find(deltafcepexp~=-1);
        delty=deltafcepexp(1,platne);
        suma=sum(delty);
        vahyw=1/suma*delty;
        mSCRpexp=max(SCRpexp);
        if MCRpexp(1,kpexp)==-1  ||  mSCRpexp==0
            MCRpexp(1,kpexp)=-1;
        else    
            MCRpexp(1,kpexp)=sum(vahyw.*SCRpexp);
        end
        meanSFpomjm=vahyw.*SFpexp;
        meanSFpomci=meanSFpomjm.*SFpexp;
        MFpexp(1,kpexp)=sum(meanSFpomci)/sum(meanSFpomjm);
        kpexp=kpexp+1;
        if kpexp>H
            kpexp=1;
        end
    end
    if uspesnychprlbin>0
        platne=find(deltafceprlbin~=-1);
        delty=deltafceprlbin(1,platne);
        suma=sum(delty);
        vahyw=1/suma*delty;
        mSCRprlbin=max(SCRprlbin);
        if MCRprlbin(1,kprlbin)==-1  || mSCRprlbin==0
            MCRprlbin(1,kprlbin)=-1;
        else    
            MCRprlbin(1,kprlbin)=sum(vahyw.*SCRprlbin);
        end
        meanSFpomjm=vahyw.*SFprlbin;
        meanSFpomci=meanSFpomjm.*SFprlbin;
        MFprlbin(1,kprlbin)=sum(meanSFpomci)/sum(meanSFpomjm);
        kprlbin=kprlbin+1;
        if kprlbin>H
            kprlbin=1;
        end
    end
    if uspesnychprlexp>0
        platne=find(deltafceprlexp~=-1);
        delty=deltafceprlexp(1,platne);
        suma=sum(delty);
        vahyw=1/suma*delty;
        mSCRprlexp=max(SCRprlexp);
        if MCRprlexp(1,kprlexp)==-1  ||  mSCRprlexp==0
            MCRprlexp(1,kprlexp)=-1;
        else    
            MCRprlexp(1,kprlexp)=sum(vahyw.*SCRprlexp);
        end
        meanSFpomjm=vahyw.*SFprlexp;
        meanSFpomci=meanSFpomjm.*SFprlexp;
        MFprlexp(1,kprlexp)=sum(meanSFpomci)/sum(meanSFpomjm);
        kprlexp=kprlexp+1;
        if kprlexp>H
            kprlexp=1;
        end
    end
    
    evals=evals+ps;
    
    ps_minule=ps;
    ps=round(((N_min-N_init)/maxiter)*evals+N_init);
    if ps<ps_minule
        P=sortrows(P,D+1);
        P=sortrows(P,D+2);
       %minimalizuju violation ale stejne tridim i podle fmin
        P=P(1:ps,:);
        max_velikost_archivu=round(ps*2.6);
        while velarchivu > max_velikost_archivu
            index_v_arch=nahvyb(velarchivu,1);
            A(index_v_arch,:)=[];
            velarchivu=velarchivu-1;
        end
    end

if (evals>=0.1*maxiter)&&(pesek1==0)
    [violmin1, indexviolmin]=min(P(:,D+2)); %minimalizuju violation
    pomocna=P(indexviolmin,:);
    [fmin1,indexfmin]=min(pomocna(:,D+1));
    bodmin1=pomocna(indexfmin,:);
    pesek1=1;
end
if (evals>=0.5*maxiter)&&(pesek2==0)
    [violmin2, indexviolmin]=min(P(:,D+2)); %minimalizuju violation
    pomocna=P(indexviolmin,:);
    [fmin2,indexfmin]=min(pomocna(:,D+1));
    bodmin2=pomocna(indexfmin,:);
    pesek2=1;
end
      
 end
[violmin, indexviolmin]=min(P(:,D+2)); %minimalizuju violation
pomocna=P(indexviolmin,:);
[fmin,indexfmin]=min(pomocna(:,D+1));
bodmin=pomocna(indexfmin,:);


c1=spocitejc(bodmin1);
c2=spocitejc(bodmin2);
c=spocitejc(bodmin);
por1=length(find(bodmin1(D+3:D+11)>0));
por2=length(find(bodmin2(D+3:D+11)>0));
por=length(find(bodmin(D+3:D+11)>0));
evals
violmin 
fmin
vystup=[evals violmin1 fmin1 c1 por1 violmin2 fmin2 c2 por2 violmin fmin c por];
