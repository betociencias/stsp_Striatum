% Synaptic Dynamics
% ODE system describing release of NT and recovery of the synapse

%function xprime = SynapticDynamics(t,x,ninf, pinf, tau_n, tau_p, pulso, q)

function xprime = SynapticDynamics(t,x, pulso)

n_inf=2;
p_inf = 0.2;
tau_n = 100;
tau_p = 10;
q = 0.2;
dx1=(1-x(1))./tau_n - pulso.*x(2).*x(1)
dx2=(p_inf-x(2))./tau_p + pulso.*p_inf.*(n_inf-x(2))
xprime = [dx1;dx2];  
%xprime = [ (ninf-x(1))./tau_n - pulso.*x(2).*x(1) ; (pinf-x(2))./tau_p + pulso.*q.*(1-x(2)) ];  
