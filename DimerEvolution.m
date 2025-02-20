% This code can be used to simulate the two-monomer (dimer) configuration
% where each site is composed of two qubits.

%First we define the maximum number of phonons

N=18;

%We will explore the single excitation manifold for 4 ions, spanned by 4
%vectors representing a single spin-up on each site.

%We construct the appropriate operators for the expectation value
%calculation

%Phonon operators:

AA = sparse(N+1,N+1);
for r=1:N
    AA(r,r+1) = sqrt(r);
end
AC = transpose(AA);
NOp = AC*AA;
Pn = sparse(N+1,N+1);
Pn(N+1,N+1)=1;

%Spin operators

Sz = sparse([-3,0,0,0;0,-3,0,0;0,0,-3,0;0,0,0,-3]);
P1 = sparse([1,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0]); %These are the projector operators into the different electronic states
P2 = sparse([0,0,0,0;0,1,0,0;0,0,0,0;0,0,0,0]);
P3 = sparse([0,0,0,0;0,0,0,0;0,0,1,0;0,0,0,0]);
P4 = sparse([0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,1]);

psiS = 1/sqrt(2)*[1,1,0,0]; %Triplet state in the donor
psiAS = 1/sqrt(2)*[-1,1,0,0]; %Singlet state in the donor
psiSAc = 1/sqrt(2)*[0,0,1,1]; %Triplet state in the acceptor
psiASAc = 1/sqrt(2)*[0,0,-1,1]; %Singlet state in the acceptor

PS = transpose(psiS)*psiS;
PAS = transpose(psiAS)*psiAS;
PSAC = transpose(psiSAc)*psiSAc;
PASAC = transpose(psiASAc)*psiASAc;

NOpT = kron(speye(4),NOp);
SzT = kron(Sz,speye(N+1));
P1T = kron(P1,speye(N+1));
P2T = kron(P2,speye(N+1));
P3T = kron(P3,speye(N+1));
P4T = kron(P4,speye(N+1));
PST=kron(PS,speye(N+1));
PAST=kron(PAS,speye(N+1));
PSACT=kron(PSAC,speye(N+1));
PASACT=kron(PASAC,speye(N+1));
PnT = kron(speye(4),Pn); %This one operator projects into the highest possible phonon number, can be used to check convergence of the phonon cutoff



%Now we set the parameters to their respective values
deltast =2*pi; %This sets all the frequencies to be angular frequencies
delta =deltast; %Phonon frequency
delz = 3*deltast; %\epsilon in the paper (single site energy)
J=0.3*deltast; %Electronic coupling
om=1*deltast; %Phonon-electronic coupling (g in the paper)
gam= 0.03*deltast; %Decay rate gamma
nbar=0.01; %Temperature of the bath (\bar{n} in the paper)
alpha = 1; %power law decay exponent  (p in the paper)

g1=om/2; g2=om/2; g3=-om/2; g4=-om/2; %We define each site phonon-spin interaction as in the paper


%Now we define all the diagonals included in our sparse matrix

diag1a=sparse((4*(N+1))^2,1);
diag1b=sparse((4*(N+1))^2,1);
diag2=sparse((4*(N+1))^2,1);
diag3=sparse((4*(N+1))^2,1);
diag4=sparse((4*(N+1))^2,1);
diag5=sparse((4*(N+1))^2,1);
diag6=sparse((4*(N+1))^2,1);
diag7=sparse((4*(N+1))^2,1);
diag8=sparse((4*(N+1))^2,1);
diag9=sparse((4*(N+1))^2,1);
diag10=sparse((4*(N+1))^2,1);
diag11=sparse((4*(N+1))^2,1);
diag12=sparse((4*(N+1))^2,1);
diag13=sparse((4*(N+1))^2,1);
diag14=sparse((4*(N+1))^2,1);
diag15=sparse((4*(N+1))^2,1);
diag16=sparse((4*(N+1))^2,1);
diag17=sparse((4*(N+1))^2,1);
diag18=sparse((4*(N+1))^2,1);
diag19=sparse((4*(N+1))^2,1);



ind1=1;ind2=1;ind3=1;ind4=1;ind5=1;ind6=1;ind7=1;ind8=1;ind9=1;ind10=1;
ind11=1;ind12=1;ind13=1;ind14=1;ind15=1;ind16=1;ind17=1;ind18=1;ind19=1;

%Now we define the generic diagonals

for q=1:4*(N+1)
    for qp =1:4*(N+1)
        if q <= N+1
            m = 1; 
            n = q-1;
        elseif q > N+1 && q<= 2*(N+1)
            m = 2;
            n = q-N-2;
        elseif q > 2*(N+1) && q<=3*(N+1)
            m= 3;
            n=q-2*N-3;
        else
            m= 4;
            n=q-3*N-4;
        end
        if qp <= N+1
            mp =1;
            np = qp-1;
        elseif qp > N+1 && qp<= 2*(N+1)
            mp =2;
            np = qp-N-2;
        elseif qp > 2*(N+1) && qp<=3*(N+1)
            mp=3;
            np=qp-2*N-3;
        else
            mp=4;
            np=qp-3*N-4;
        end
        Aan = An(n);
        Bbn = Bn(n,N);
        Aanp = An(np);
        Bbnp = Bn(np,N);

        %The diagonal terms

        diag1a(ind1,1) = -1i*(delta*n-delta*np) -gam*(nbar+1)*(n+np)/2 - gam*nbar*(n+np+2)/2;
        diag1b(ind1,1) = -1i*(kronDel(m,1)+kronDel(m,2)-kronDel(m,3)-kronDel(m,4) - kronDel(mp,1)-kronDel(mp,2)+kronDel(mp,3)+kronDel(mp,4));
        ind1=ind1+1;

        %The spin-phonon terms

        diag2(ind2,1)= -1i*Aan*(g1*kronDel(m,1)+g2*kronDel(m,2)+g3*kronDel(m,3)+g4*kronDel(m,4));
        ind2=ind2+1;

        diag3(ind3,1)= -1i*Bbn*(g1*kronDel(m,1)+g2*kronDel(m,2)+g3*kronDel(m,3)+g4*kronDel(m,4));
        ind3=ind3+1;

        diag4(ind4,1)= 1i*Aanp*(g1*kronDel(mp,1)+g2*kronDel(mp,2)+g3*kronDel(mp,3)+g4*kronDel(mp,4));
        ind4=ind4+1;

        diag5(ind5,1)= 1i*Bbnp*(g1*kronDel(mp,1)+g2*kronDel(mp,2)+g3*kronDel(mp,3)+g4*kronDel(mp,4));
        ind5=ind5+1;
        
       %The spin-spin terms

        diag6(ind6,1)= -1i*J*(kronDel(m,2)+(1/3^alpha)*kronDel(m,3)+kronDel(m,4));
        ind6=ind6+1;

        diag7(ind7,1)= -1i*J*(kronDel(m,1)+(1/3^alpha)*kronDel(m,2)+kronDel(m,3));
        ind7=ind7+1;

        diag8(ind8,1)= -1i*(J/(4^alpha))*(kronDel(m,3)+kronDel(m,4)); 
        ind8=ind8+1;

        diag9(ind9,1)= -1i*(J/(4^alpha))*(kronDel(m,1)+kronDel(m,2)); 
        ind9=ind9+1;
        
        diag10(ind10,1) = -1i*(J/(5^alpha))*kronDel(m,4); 
        ind10 = ind10+1;

        diag11(ind11,1) = -1i*(J/(5^alpha))*kronDel(m,1); 
        ind11 = ind11+1;

        diag12(ind12,1) = 1i*J*(kronDel(mp,2)+(1/3^alpha)*kronDel(mp,3)+ kronDel(mp,4));
        ind12 = ind12+1;

        diag13(ind13,1) = 1i*J*(kronDel(mp,1)+(1/3^alpha)*kronDel(mp,2)+kronDel(mp,3));
        ind13 = ind13+1;

        diag14(ind14,1) = 1i*(J/(4^alpha))*(kronDel(mp,3)+kronDel(mp,4));
        ind14 = ind14+1;

        diag15(ind15,1) = 1i*(J/(4^alpha))*(kronDel(mp,1)+kronDel(mp,2));
        ind15 = ind15+1;

        diag16(ind16,1) = 1i*(J/(5^alpha))*kronDel(mp,4);
        ind16 = ind16+1;

        diag17(ind17,1) = 1i*(J/(5^alpha))*kronDel(mp,1); 
        ind17 = ind17+1;
        

        %Dissipation terms
        
        diag18(ind18,1) = gam*(nbar+1)*Bbn*Bbnp;
        ind18=ind18+1;

        diag19(ind19,1) = gam*(nbar)*Aan*Aanp;
        ind19 = ind19+1;


    end
end

TransferC = zeros(1,3);
indT = 1;

%Now we define the matrix of coefficients

SSDM = spdiags([diag1a+delz/2*diag1b diag4 diag5 diag2 diag3 diag12 diag13 diag14 diag15 diag16 diag17 diag6 diag7 diag8 diag9 diag10 diag11 diag18 diag19],[0 1 -1 4*(N+1) -4*(N+1) N+1 -N-1 2*N+2 -2*N-2 3*N+3 -3*N-3 4*(N+1)*(N+1) 4*(N+1)*(-N-1) 4*(N+1)*(2*N+2) 4*(N+1)*(-2*N-2) 4*(N+1)*(3*N+3) 4*(N+1)*(-3*N-3) -4*N-5 4*N+5],(4*(N+1))^2,(4*(N+1))^2);

%Note that in the line above we multiply here by delz/2 such that the user
%can create a for loop with different values of delz (\epsilon in the
%paper) in case they want to create figures similar to fig.2 (a) and (c) in
%the paper

SSDM = transpose(SSDM);

y0 = reshape(thermalState(N,nbar,om),(4*(N+1))^2,1); %The initial state is defined with the function thermalState described below

tspan = linspace(0,100,300); %Time for the simulation (immediately in units of t \omega / 2pi)

[t,y] = ode45(@(t,y) ODEDiss(t,y,SSDM),tspan,y0); %Stores the density matrix for each step

results = zeros(50,14); %Array to store observables

for jj=1:size(y,1)

SS = reshape(y(jj,:),4*(N+1),4*(N+1));

results(jj,1) = t(jj);
results(jj,2) = trace(NOpT*SS);
results(jj,3) = trace(P1T*SS);
results(jj,4) = trace(P2T*SS);
results(jj,5) = trace(P3T*SS);
results(jj,6) = trace(P4T*SS);
results(jj,7) = trace(PnT*SS);
results(jj,8) = trace(SS*SS);
results(jj,12) = trace(PST*SS);
results(jj,13) = trace(PSACT*SS);
results(jj,14) = trace(PASACT*SS);
end

%%Finding the transfer rate on Eq. 14 on the paper

Sum1 = trapz(tspan,results(:,4)+results(:,3)); %numerator of the expression
Sum2 = trapz(tspan,transpose(tspan).*(results(:,4)+results(:,3)));% denominator of the expression
kappaT = Sum1/Sum2; %k_T in the paper


%%
function dy = ODEDiss(t,y,M)
dy = M*y;
end

function x = An(n)
    x = sqrt(n);
end

function x = Bn(n,N)
if n >= N
    x = 0;
else
    x = sqrt(n+1);
end
end

function d=kronDel(j,v)
d=j==v;
end

function rho=thermalState(N,nbar,g)
%This functions defines a given initial state with the phonon state being a
%thermal state determined by its temperature nbar and displacement g

AA = sparse(N+1,N+1);
for r=1:N
    AA(r,r+1) = sqrt(r);
end
AC = transpose(AA);
NOp = AC*AA;

%We generate the displaced vacuum
boson=sparse(N+1,1);
alpha = g/2;
prefac = exp(-(alpha)^2/2);
for jj=0:N
    boson(jj+1) = prefac*(alpha^(jj)/(sqrt(factorial(jj))));
end
boson=boson/(norm(boson));

if nbar ~= 0
    kbt = 1/(log(1+1/nbar));
    %spin state
    state1 = [1;1;0;0]; %Here we define the electronic state to be the triplet in the donor but this can be changed
    state1=state1/(norm(state1));
    rhospin = sparse(state1*ctranspose(state1));
    %Bosonic dm
    rhob = sparse(N+1,N+1);
    for n=1:N+1
        newbos = 1/(sqrt(factorial(n-1)))*((AC+alpha*speye(N+1))^(n-1))*boson;
        newbos = newbos/norm(newbos);
        coeff = 1/(nbar+1)*(nbar/(nbar+1))^(n-1);
        rhob = rhob+coeff*newbos*ctranspose(newbos);
    end
    trace(rhob);
    rho = kron(rhospin,rhob);
    trace(rhob*NOp);
    nbar+alpha^2;
else
    state1 = [1;1;0;0]; %Here we define the electronic state to be the triplet in the donor but this can be changed
    state1=state1/(norm(state1));
    statefull = kron(state1,boson);
    rho = statefull*ctranspose(statefull);
end


end