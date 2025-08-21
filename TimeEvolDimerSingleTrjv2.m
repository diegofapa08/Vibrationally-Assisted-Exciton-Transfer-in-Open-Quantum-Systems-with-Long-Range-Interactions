% Now we will build a kind of artificial Hamiltonian just to understand
% How to define the different transport regimes

%First we define the maximum number of phonons

N=18;
Nelec = 4;
dim = (N+1)*Nelec;
%We will explore the single excitation manifold for 4 ions, spanned by 4
%vectors representing a single spin-up on each site.

%We construct the appropriate operators for the expectation value
%calculation

AA = sparse(N+1,N+1);
for r=1:N
    AA(r,r+1) = sqrt(r);
end
AC = transpose(AA);
NOp = AC*AA;
Pn = sparse(N+1,N+1);
Pn(N+1,N+1)=1;

alpha =1;
Jmat = sparse(eCoupling(Nelec,2,2,alpha)); %Jij matrix of electronic couplings
gmat = sparse(ebCoupling(Nelec,2)); %coefficients of g(a^dag + a) for each electronic state
epsmat = sparse(eShift(Nelec,2)); %shifts for each state

P1 = sparse([1,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0]);
P2 = sparse([0,0,0,0;0,1,0,0;0,0,0,0;0,0,0,0]);
P3 = sparse([0,0,0,0;0,0,0,0;0,0,1,0;0,0,0,0]);
P4 = sparse([0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,1]);

NOpT = kron(speye(4),NOp);
P1T = kron(P1,speye(N+1));
P2T = kron(P2,speye(N+1));
P3T = kron(P3,speye(N+1));
P4T = kron(P4,speye(N+1));
PnT = kron(speye(4),Pn);

%Hamiltonian parameters 

delta = 2*pi;
J=0.3*delta;
g = delta;
gam = 0.039552*delta;
nbar =0.01;
epsilon = 3*delta;


Hm = delta*NOpT + g/2*kron(gmat,AA+AC) + epsilon/2*kron(epsmat,speye(N+1)) + J*kron(Jmat,speye(N+1));

JopM = kron(speye(Nelec),AA);
JopP = JopM';

LOP = -1i*(kron(speye(dim),Hm)-kron(transpose(Hm),speye(dim)));
LOP = LOP + gam*(nbar+1)*(kron(JopM,conj(JopM))-kron(speye(dim),JopM'*JopM)/2 -kron(transpose(JopM'*JopM),speye(dim))/2);
LOP = LOP + gam*(nbar)*(kron(JopP,conj(JopP))-kron(speye(dim),JopP'*JopP)/2 -kron(transpose(JopP'*JopP),speye(dim))/2);

%Initial state and simulation details
y0 = reshape(thermalState(N,nbar,g,delta),(Nelec*(N+1))^2,1);
totaltime = 100; %In the units indicated in the x-axis of the figures of the paper, number of cycles.
Nsteps = 300;
deltat = totaltime/Nsteps;
t=-deltat;
results = zeros(50,14);
statef=y0;
for jj=1:Nsteps+1
    t = t+deltat;
    SS = reshape(statef,dim,dim);
    results(jj,1) = t;
    results(jj,2) = trace(NOpT*SS);
    results(jj,3) = trace(P1T*SS);
    results(jj,4) = trace(P2T*SS);
    results(jj,5) = trace(P3T*SS);
    results(jj,6) = trace(P4T*SS);
    results(jj,7) = trace(PnT*SS);
    results(jj,8) = trace(SS*SS);
    statef = expv(deltat,LOP,statef);
    jj
end
%%
tspan = results(:,1);
Sum1 = trapz(tspan,results(:,4)+results(:,3));
Sum2 = trapz(tspan,(tspan).*(results(:,4)+results(:,3)));
kappaT = Sum1/Sum2;
kappaT = (kappaT);

%%
plot(results(:,1),results(:,3)+results(:,4))
hold on
plot(results(:,1),results(:,5)+results(:,6))
plot(results(:,1),results(:,3)+results(:,5)+results(:,4)+results(:,6));
xlim([0,50])
hold off

%%
function dy = ODEDiss(t,y,M)
dy = M*y;
end

function Jij = eCoupling(N,d,m,al)
%Power law decaying coupling: N is the total number of electronic sites, d
%is the distance between monomers, m is the number of qubits per
%monomer, and al is the power law exponent
Jij = zeros(N,N);

states_per_monomer = m;
subunit_labels = ceil((1:N) / states_per_monomer);  % vector of subunit indices


% Fill the matrix with modified distances
for i = 1:N
    for j = 1:N
        if i ~= j
            base_dist = abs(i - j);
            subunit_dist = abs(subunit_labels(i) - subunit_labels(j));
            total_dist = base_dist + d * subunit_dist;
            Jij(i,j) = 1 / total_dist^al;
        end
    end
end


end

function gmat = ebCoupling(N,m)
%Electronic-Boson coupling: N is the total number of electronic sites, m is the number of qubits per
%monomer.
gmat= zeros(N,N);

states_per_monomer = m;
subunit_labels = ceil((1:N) / states_per_monomer);  % vector of subunit indices


% Fill the matrix with modified distances
for i = 1:N
    coeff = subunit_labels(i);
    gmat(i,i) = (-1)^(1+coeff);
end


end

function epsmat = eShift(N,m)
%Electronic-Boson coupling: N is the total number of electronic sites, m is the number of qubits per
%monomer.
epsmat= zeros(N,N);

states_per_monomer = m;
subunit_labels = ceil((1:N) / states_per_monomer);  % vector of subunit indices


% Fill the matrix with modified distances
for i = 1:N
    coeff = subunit_labels(i);
    epsmat(i,i) = 3-2*coeff;
end


end


function rho=thermalState(N,nbar,g,om)

AA = sparse(N+1,N+1);
for r=1:N
    AA(r,r+1) = sqrt(r);
end
AC = transpose(AA);
NOp = AC*AA;

%We generate the displaced vacuum
boson=sparse(N+1,1);
alpha = g/(2*om);
prefac = exp(-(alpha)^2/2);
for jj=0:N
    boson(jj+1) = prefac*(alpha^(jj)/(sqrt(factorial(jj))));
end
boson=boson/(norm(boson));

if nbar ~= 0
    kbt = 1/(log(1+1/nbar));
    %spin state
    state1 = [1;1;0;0];
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
    state1 = [1;1;0;0];
    state1=state1/(norm(state1));
    statefull = kron(state1,boson);
    rho = statefull*ctranspose(statefull);
end


end