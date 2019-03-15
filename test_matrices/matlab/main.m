%%
%{
The goal of this script is to investigate the convergence of the matlab implementations
of BiCGStab and BiCGStab(l). By default l=2, no preconditioner is used and
an initial guess is the zero-vector. 
%}

%% Case 1
%{
Discretization of the differential equation:
-u_xx-u_yy+1000*(x*u_x+y*u_y)+10u=f
where the solution is defined such that the solution vector of u is all
ones. 
This differential equation is a test problem in the paper:
Sleijpen, G. L., Van der Vorst, H. A., & Fokkema, D. R. (1994). 
BiCGstab(l) and other hybrid Bi-CG methods. Numerical Algorithms, 7(1), 75-109. doi:10.1007/bf02141261


The number of iterations is slightly weird, matlab's BiCGStab performs
"half iterations", BiCGStab(l) performs "quarter iterations". Each fractional
iteration yields a residual. 
By default matlab's BiCGStab(l) uses l=2.
%}

clear all
close all
format compact

Lx=1;
Nx=64;
beta=1000;

Ny=Nx;
Ly=Lx;

dx=Lx./(Nx);
dx2=dx.*dx;
dy=Ly./(Ny);
dy2=dy.*dy;

N=Nx;

A=sparse(Nx,Ny,5.*Nx);
b=ones(Nx.*Ny,1).*10;

%{
Matrix assembly
%}
for i=0:N-1 %Loop over grid
    for j=0:N-1
        current_point=j*N+i;
        coeff_sum=0;
        
        A(current_point+1,current_point+1)=(2./dx2+2./dy2+10);
        coeff_sum=coeff_sum+A(current_point+1,current_point+1);
        %north=[i-1,j];
        k=j*N+(i-1);
        if i>0
            A(current_point+1,k+1)=(-1./dy2-beta.*(dy./2+i.*dy)./(2.*dy));
            coeff_sum=coeff_sum+A(current_point+1,k+1);
        end
        %east=[i,j+1];
        k=(j+1)*N+i;
        if j<N-1
            A(current_point+1,k+1)=(-1./dx2+beta.*(dx./2+j.*dx)./(2.*dx));
            coeff_sum=coeff_sum+A(current_point+1,k+1);
        end
        %south=[i+1,j];
        k=j*N+(i+1);
        if i<N-1
            A(current_point+1,k+1)=(-1./dy2+beta.*(dy./2+i.*dy)./(2.*dy));
            coeff_sum=coeff_sum+A(current_point+1,k+1);
        end
        %west=[i,j-1];
        k=(j-1)*N+i;
        if j>0
            A(current_point+1,k+1)=(-1./dx2-beta.*(dx./2+j.*dx)./(2.*dx));
            coeff_sum=coeff_sum+A(current_point+1,k+1);
        end
        b(current_point+1)=coeff_sum; %Ensure that the exact solution vector will be all ones
    end
end

tol=1e-12;
maxit=2000;

%{
Compute with BiCGStab & BiCGStab(2)
%}
[x_stab,flag_stab,relres_stab,iter_stab,resvec_stab] = bicgstab(A,b,tol,maxit);
[x_stabL,flag_stabL,relres_stabL,iter_stabL,resvec_stabL] = bicgstabl(A,b,tol,maxit);

%{
flags: 
0 converged to the desired tolerance within the set number of iterations
1 reached maximum number of iterations, but did not converge
2 precondtioner was ill-conditionerd
3 stagnated (two consecutive iterates were the same)
4 one of the scalar quantities calculated became too small or too large to
continue computing
%}
flag_stab
flag_stabL

%{
Make plots
%}
figure(2)
semilogy([0.5:0.5:length(resvec_stab)/2],resvec_stab,'displayname','Residual BiCGStab')
hold on
semilogy([0.25:0.25:length(resvec_stabL)/4],resvec_stabL,'displayname','Residual BiCGStab(2)')
xlabel('Iteration')
ylabel('||Ax-b||_2')
legend('-dynamiclegend','location','southwest')
axis tight
set(gca,'linewidth',1.5)

%{
Save plot to file for use in gnuplot
%}
M=[transpose([0.5:0.5:length(resvec_stab)/2]) resvec_stab];
dlmwrite('bicgstab_sleijpen.txt',M)
M=[transpose([0.25:0.25:length(resvec_stabL)/4]) resvec_stabL];
dlmwrite('bicgstab2_sleijpen.txt',M)

eigs(A,10)
condest(A)

%{
Export A,b,x in matrix market format
%}
mmwrite('sleijpen_example_A.txt',A,'Discretization matrix from G.L.G. Sleijpen BiCGstab(l) paper.')
mmwrite('sleijpen_example_b.txt',b,'Discretization matrix from G.L.G. Sleijpen BiCGstab(l) paper.')
x=A\b;
mmwrite('sleijpen_example_x.txt',x,'Discretization matrix from G.L.G. Sleijpen BiCGstab(l) paper.')

%% Case 2
%{
Construct a simple amtrix for which bicgstab does not converge, but
bicgstabl does. The original idea was to derive some analytical results. 

One can show 2 things:
1. This is a so called "tridiagonal Toeplitz matrix". The eigenvalues of
    the matrix A are b+2sqrt(-a^2)cos(k*pi/(N+1)) where k is an integer [1,...N]. 
2. The parameter omega in BiCGStab is suspected to scale ~c1/(c2+a^2) for
    large a and some constants c1,c2. However the algorithm did not seem to 
    stop because of omega in tests. 
%}

N=200;
N2=N*N; %A is an N^2xN^2 matrix.
a=10;
b=1;
%{
For sufficiently large a the matrix becomes strongly asymmetric.
Large b (diagonal elements) can sort-of alleviate problems. 
The general trend seems to be:
-Larger a -> bicgstab can have problems
-Extremely large a, both have problems
-Larger b -> convergence much nicer for both

Known case where bicgstab diverges, but bicgstab(l) converges:
N=100,200
a=10
b=1
%}

%{
Matrix assembly
%}
vec=(ones(N2,1))*b;
A = spdiags(vec(:),0,N2,N2);
vec=(ones(N2,1))*a;
A = spdiags(vec,1,A);
vec=(ones(N2,1))*(-a);
A = spdiags(vec,-1,A);
b=ones(N2,1);

%{
Compute with BiCGStab & BiCGStab(2)
%}
[x_stab,flag_stab,relres_stab,iter_stab,resvec_stab] = bicgstab(A,b,tol,maxit);
[x_stabL,flag_stabL,relres_stabL,iter_stabL,resvec_stabL] = bicgstabl(A,b,tol,maxit);

%{
flags: 
0 converged to the desired tolerance within the set number of iterations
1 reached maximum number of iterations, but did not converge
2 precondtioner was ill-conditionerd
3 stagnated (two consecutive iterates were the same)
4 one of the scalar quantities calculated became too small or too large to
continue computing
%}
flag_stab
flag_stabL

%{
Make plots
%}
figure(4)
semilogy([0.5:0.5:length(resvec_stab)/2],abs(resvec_stab),'displayname','Residual BiCGStab')
hold on
semilogy([0.25:0.25:length(resvec_stabL)/4],abs(resvec_stabL),'displayname','Residual BiCGStab(2)')
xlabel('Iteration')
ylabel('||Ax-b||_2')
legend('-dynamiclegend','location','southeast')
title([''])
axis tight

eigs(A,10)
condest(A)

%{
Save plot to file for use in gnuplot
%}
M=[transpose([0.5:0.5:length(resvec_stab)/2]) resvec_stab];
dlmwrite(['bicgstab_tridiagonal_a_',num2str(a),'.txt'],M)
M=[transpose([0.25:0.25:length(resvec_stabL)/4]) resvec_stabL];
dlmwrite(['bicgstab2_tridiagonal_a_',num2str(a),'.txt'],M)

%{
Export A,b,x in matrix market format
%}
mmwrite('toeplitz_A.txt',A,'Tridiagonal Toeplitz test matrix.')
mmwrite('toeplitz_b.txt',b,'Tridiagonal Toeplitz test matrix.')
x=A\b;
mmwrite('toeplitz_x.txt',x,'Tridiagonal Toeplitz test matrix.')



