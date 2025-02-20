% This code can be used to simulate the three-monomer configuration
% where each site is composed of two qubits.

%First we define the maximum number of phonons

N=18;

%We will explore the single excitation manifold for 6 ions, spanned by 6
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

Sz = sparse([-5,0,0,0,0,0;0,-5,0,0,0,0;0,0,-5,0,0,0;0,0,0,-5,0,0;0,0,0,0,-5,0;0,0,0,0,0,-5]);
P1 = sparse([1,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0]); %These are projectors into each electronic state
P2 = sparse([0,0,0,0,0,0;0,1,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0]);
P3 = sparse([0,0,0,0,0,0;0,0,0,0,0,0;0,0,1,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0]);
P4 = sparse([0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,1,0,0;0,0,0,0,0,0;0,0,0,0,0,0]);
P5= sparse([0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,1,0;0,0,0,0,0,0]);
P6=sparse([0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,1]);

NOpT = kron(speye(6),NOp);
SzT = kron(Sz,speye(N+1));
P1T = kron(P1,speye(N+1));
P2T = kron(P2,speye(N+1));
P3T = kron(P3,speye(N+1));
P4T = kron(P4,speye(N+1));
P5T = kron(P5,speye(N+1));
P6T = kron(P6,speye(N+1));
PnT = kron(speye(6),Pn); %This one operator projects into the highest possible phonon number, can be used to check convergence of the phonon cutoff



%Now we set the parameters to their respective values
deltast =2*pi; %This sets all the frequencies to be angular frequencies
delta =deltast; %Phonon frequency
delz = 3*deltast; %\epsilon in the paper (single site energy)
J=0.3*deltast; %Electronic coupling
om=1*deltast; %Phonon-electronic coupling (g in the paper)
gam= 0.03*deltast; %Decay rate gamma
nbar=0.01; %Temperature of the bath (\bar{n} in the paper)
alpha = 1; %power law decay exponent  (p in the paper)

g1=om/2;g2=om/2;g3=-om/2;g4=-om/2;g5=om/2;g6=om/2; %We define each site phonon-spin interaction as in the paper


for q=1:6*(N+1)
    for qp =1:6*(N+1)
        if q <= N+1
            m = 1; 
            n = q-1;
        elseif q > N+1 && q<= 2*(N+1)
            m = 2;
            n = q-N-2;
        elseif q > 2*(N+1) && q<=3*(N+1)
            m= 3;
            n=q-2*N-3;
        elseif q > 3*(N+1) && q<=4*(N+1)
            m= 4;
            n=q-3*N-4;
        elseif q > 4*(N+1) && q<=5*(N+1)
            m= 5;
            n=q-4*N-5;
        else
            m= 6;
            n=q-5*N-6;
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
        elseif qp > 3*(N+1) && qp<=4*(N+1)
            mp=4;
            np=qp-3*N-4;
        elseif qp > 4*(N+1) && qp<=5*(N+1)
            mp=5;
            np=qp-4*N-5;
        else
            mp=6;
            np=qp-5*N-6;
        end
        Aan = An(n);
        Bbn = Bn(n,N);
        Aanp = An(np);
        Bbnp = Bn(np,N);

        %The diagonal terms

        diag1a(ind1,1) = -1i*(delta*n-delta*np) -gam*(nbar+1)*(n+np)/2 - gam*nbar*(n+np+2)/2;
        diag1b(ind1,1) = -1i*(kronDel(m,1)+kronDel(m,2)-kronDel(m,3)-kronDel(m,4)-3*kronDel(m,5)-3*kronDel(m,6)-kronDel(mp,1)-kronDel(mp,2)+kronDel(mp,3)+kronDel(mp,4)+3*kronDel(mp,5)+3*kronDel(mp,6));
        ind1=ind1+1;

        %The spin-phonon terms

        diag2(ind2,1)= -1i*Aan*(g1*kronDel(m,1)+g2*kronDel(m,2)+g3*kronDel(m,3)+g4*kronDel(m,4)+g5*kronDel(m,5)+g6*kronDel(m,6));
        ind2=ind2+1;

        diag3(ind3,1)= -1i*Bbn*(g1*kronDel(m,1)+g2*kronDel(m,2)+g3*kronDel(m,3)+g4*kronDel(m,4)+g5*kronDel(m,5)+g6*kronDel(m,6));
        ind3=ind3+1;

        diag4(ind4,1)= 1i*Aanp*(g1*kronDel(mp,1)+g2*kronDel(mp,2)+g3*kronDel(mp,3)+g4*kronDel(mp,4)+g5*kronDel(mp,5)+g6*kronDel(mp,6));
        ind4=ind4+1;

        diag5(ind5,1)= 1i*Bbnp*(g1*kronDel(mp,1)+g2*kronDel(mp,2)+g3*kronDel(mp,3)+g4*kronDel(mp,4)+g5*kronDel(mp,5)+g6*kronDel(mp,6));
        ind5=ind5+1;
        
       %The spin-spin terms

        diag6(ind6,1)= -1i*J*(kronDel(m,2)+(1/3^alpha)*kronDel(m,3)+kronDel(m,4)+(1/3^alpha)*kronDel(m,5)+kronDel(m,6));
        ind6=ind6+1;

        diag7(ind7,1)= -1i*J*(kronDel(m,1)+(1/3^alpha)*kronDel(m,2)+kronDel(m,3)+(1/3^alpha)*kronDel(m,4)+kronDel(m,5));
        ind7=ind7+1;

        diag8(ind8,1)= -1i*(J/(4^alpha))*(kronDel(m,3)+kronDel(m,4)+kronDel(m,5)+kronDel(m,6)); 
        ind8=ind8+1;

        diag9(ind9,1)= -1i*(J/(4^alpha))*(kronDel(m,1)+kronDel(m,2)+kronDel(m,3)+kronDel(m,4)); 
        ind9=ind9+1;
        
        diag10(ind10,1) = -1i*(J/(5^alpha))*(kronDel(m,4)+kronDel(m,6))-1i*(J/(7^alpha))*(kronDel(m,5)); 
        ind10 = ind10+1;

        diag11(ind11,1) = -1i*(J/(5^alpha))*(kronDel(m,1)+kronDel(m,3))-1i*(J/(7^alpha))*(kronDel(m,2)); 
        ind11 = ind11+1;

        diag12(ind12,1) = -1i*(J/(8^alpha))*(kronDel(m,6)+kronDel(m,5)); 
        ind12 = ind12+1;

        diag13(ind13,1) = -1i*(J/(8^alpha))*(kronDel(m,2)+kronDel(m,1)); 
        ind13 = ind13+1;

        diag14(ind14,1) = -1i*(J/(9^alpha))*(kronDel(m,6));  
        ind14 = ind14+1;

        diag15(ind15,1) = -1i*(J/(9^alpha))*(kronDel(m,1)); %%
        ind15 = ind15+1;

        diag16(ind16,1)= 1i*J*(kronDel(mp,2)+(1/3^alpha)*kronDel(mp,3)+kronDel(mp,4)+(1/3^alpha)*kronDel(mp,5)+kronDel(mp,6));
        ind16=ind16+1;

        diag17(ind17,1)= 1i*J*(kronDel(mp,1)+(1/3^alpha)*kronDel(mp,2)+kronDel(mp,3)+(1/3^alpha)*kronDel(mp,4)+kronDel(mp,5));
        ind17=ind17+1;

        diag18(ind18,1)= 1i*(J/(4^alpha))*(kronDel(mp,3)+kronDel(mp,4)+kronDel(mp,5)+kronDel(mp,6)); 
        ind18=ind18+1;

        diag19(ind19,1)= 1i*(J/(4^alpha))*(kronDel(mp,1)+kronDel(mp,2)+kronDel(mp,3)+kronDel(mp,4)); 
        ind19=ind19+1;
        
        diag20(ind20,1) = 1i*(J/(5^alpha))*(kronDel(mp,4)+kronDel(mp,6))+1i*(J/(7^alpha))*(kronDel(mp,5)); 
        ind20 = ind20+1;

        diag21(ind21,1) = 1i*(J/(5^alpha))*(kronDel(mp,1)+kronDel(mp,3))+1i*(J/(7^alpha))*(kronDel(mp,2)); 
        ind21 = ind21+1;

        diag22(ind22,1) = 1i*(J/(8^alpha))*(kronDel(mp,6)+kronDel(mp,5)); 
        ind22 = ind22+1;

        diag23(ind23,1) = 1i*(J/(8^alpha))*(kronDel(mp,2)+kronDel(mp,1)); 
        ind23 = ind23+1;

        diag24(ind24,1) = 1i*(J/(9^alpha))*(kronDel(mp,6));  
        ind24 = ind24+1;

        diag25(ind25,1) = 1i*(J/(9^alpha))*(kronDel(mp,1)); %%
        ind25 = ind25+1;

        %Dissipation terms
        
        diag26(ind26,1) = gam*(nbar+1)*Bbn*Bbnp;
        ind26=ind26+1;

        diag27(ind27,1) = gam*(nbar)*Aan*Aanp;
        ind27 = ind27+1;


    end
end


%Now we define the matrix of coefficients

SSDM = spdiags([diag1a+delz/2*diag1b diag4 diag5 diag2 diag3 diag16 diag17 diag18 diag19 diag20 diag21 diag22 diag23 diag24 diag25 diag6 diag7 diag8 diag9 diag10 diag11 diag12 diag13 diag14 diag15 diag26 diag27],[0 1 -1 6*(N+1) -6*(N+1) N+1 -N-1 2*N+2 -2*N-2 3*N+3 -3*N-3 4*N+4 -4*N-4 +5*N+5 -5*N-5 6*(N+1)*(N+1) 6*(N+1)*(-N-1) 6*(N+1)*(2*N+2) 6*(N+1)*(-2*N-2) 6*(N+1)*(3*N+3) 6*(N+1)*(-3*N-3) 6*(N+1)*(4*N+4) 6*(N+1)*(-4*N-4) 6*(N+1)*(5*N+5) 6*(N+1)*(-5*N-5) -6*N-7 6*N+7],(6*(N+1))^2,(6*(N+1))^2);

%Note that in the line above we multiply here by delz/2 such that the user
%can create a for loop with different values of delz (\epsilon in the
%paper) in case they want to create figures similar to fig.2 (a) and (c) in
%the paper

SSDM = transpose(SSDM);

y0 = reshape(thermalStateTrim(N,nbar,om),(6*(N+1))^2,1); %The initial state is defined with the function thermalState described below

tspan = linspace(0,10,10); %Time for the simulation (immediately in units of t \omega / 2pi)

[t,y] = ode45(@(t,y) ODEDiss(t,y,SSDM),tspan,y0); %Stores the density matrix for each step

results = zeros(10,14); %Array to store observables

for jj=1:size(y,1)

SS = reshape(y(jj,:),6*(N+1),6*(N+1));

results(jj,1) = t(jj);
results(jj,2) = trace(NOpT*SS);
results(jj,3) = trace(P1T*SS);
results(jj,4) = trace(P2T*SS);
results(jj,5) = trace(P3T*SS);
results(jj,6) = trace(P4T*SS);
results(jj,7) = trace(P5T*SS);
results(jj,8) = trace(P6T*SS);
results(jj,9) = trace(PnT*SS);
results(jj,10) = trace(SS*SS);

end

%%Finding the transfer rate on Eq. 14 on the paper
%%
Sum1=0;
Sum2=0;
Sum1 = trapz(tspan,1-(results(:,7)+results(:,8))); %Numerator of the expression
Sum2 = trapz(tspan,transpose(tspan(1:end)).*(1-(results(:,7)+results(:,8)))); %Denominator of the expression
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

function rho=thermalStateTrim(N,nbar,g)
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
    state1 = [1;1;0;0;0;0]; %Here we define the electronic state to be the triplet in the donor but this can be changed
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
    state1 = [1;1;0;0;0;0]; %Here we define the electronic state to be the triplet in the donor but this can be changed
    state1=state1/(norm(state1));
    statefull = kron(state1,boson);
    rho = statefull*ctranspose(statefull);
end


end