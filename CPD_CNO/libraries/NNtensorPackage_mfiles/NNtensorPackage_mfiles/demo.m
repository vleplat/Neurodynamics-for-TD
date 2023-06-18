s1 = input('First tensor dimension: \n'); if s1<2,  s1=2; elseif s1>10, s1=10;end;
s2 = input('Second tensor dimension: \n'); if s2<2, s2=2; elseif s2>8,  s2=8; end;
s3 = input('Third tensor dimension: \n');if s3<2,   s3=2; elseif s3>6,  s3=6 ;end;

nbFactors = input('Number of factors: \n'); 
if nbFactors < 2, nbFactors = 2; 
elseif nbFactors > min(s1,(min(s2, s3))), nbFactors = min(s1,(min(s2, s3))); end;

%s1 = 2;%10;
%s2 = 3;%9;
%s3 = 4;%6;
%nbFactors = 2;%3;

rand('twister',sum(100*clock))

disp('Let construct a random tensor, and then estimate the loadings')
%initialize final positive random loadings
Af = rand(s1, nbFactors);
Bf = rand(s2, nbFactors);
Cf = rand(s3, nbFactors);

%Construct a positive 3rd tensor with these loadings
T = Af * khatriRao(Cf, Bf)';
T = reshape(T, s1, s2, s3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now, apply some algorithms to recover Af, Bf, Cf from T.

AInit = rand(size(Af));
BInit = rand(size(Bf));
CInit = rand(size(Cf));

%First, the conjugate gradient
disp('First, the conjugate gradient. Press a key to continue')
pause
%Used with a random initialization
%[Acg Bcg C errorcg] = cgP(T, nbFactors);
[Acg Bcg C errorcg] = cgP(T, nbFactors, [], AInit, BInit, CInit);

%Next, the gradient
disp('Next, the gradient. Press a key to continue')
pause
%Used with a random initialization
%[Ag Bg Cg errorg] = gradP(T, nbFactors);
[Ag Bg Cg errorg] = gradP(T, nbFactors, [], AInit, BInit, CInit);

%Finally, the bfgs
disp('Finally, the BFGS. Press a key to continue')
pause
%Used with a random initialization
%[Ab Bb Cb errorb] = bfgsP(T, nbFactors);
[Ab Bb Cb errorb] = bfgsP(T, nbFactors, [], AInit, BInit, CInit);

%display the results
disp('Let see the results.')
figure;
hold on
plot(10 * log10(errorcg))
plot(10 * log10(errorg), 'r')
plot(10 * log10(errorb), 'g')
legend('CG', 'Gradient', 'BFGS')
title('Reconstruction error')
xlabel('Iterations')
ylabel('Reconstruction error')

%We can add options when lauching the algorithms.
%See createOptions for more informations
%For example:
disp(' ')
disp(' ')
disp('Example with some options... (see createOptions)')
disp(' ')
[A B C error] = cgP(T, nbFactors, createOptions(1, 0.0001, 0, 100));

disp(' ')
disp(' ')
disp('Other example with a fixed initialization:')
disp(' ')

%Note that the initial loadings don't have to be positive (and it is not a
%problem). But the final ones estimated WILL be positive.

A = rand(size(Af));
B = rand(size(Bf));
C = rand(size(Cf));

[A B C error] = cgP(T, nbFactors, createOptions(1, 0.0001, 0, 100), A, B, C);
