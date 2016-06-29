% Generate linear regression data set with outliers
nInstances = 400;
nVars = 1;
[X,y] = makeData('regressionOutliers',nInstances,nVars);

% Least squares solution
wLS = X\y;

options = [];
options.display = 'none';
%options.maxFunEvals = maxFunEvals;
options.Method = 'lbfgs';

% Huber loss
changePoint = .2;
fprintf('Training Huber robust regression model...\n');
wHuber = minFunc(@HuberLoss,wLS,options,X,y,changePoint);

% Student T loss
lambda = 1;
dof = 2;
funObj = @(params)studentLoss(X,y,params(1:nVars),params(nVars+1),params(end));
fprintf('Training student T robust regression model...\n');
params = minFunc(funObj,[wLS;lambda;dof],options);
wT = params(1:nVars);
lambda = params(nVars+1);
dof = params(end);

% Plot results
figure;hold on
plot(X,y,'.');
xl = xlim;
h1=plot(xl,xl*wLS,'r');
h2=plot(xl,xl*wHuber,'g');
h3=plot(xl,xl*wT,'k--');
set(h1,'LineWidth',3);
set(h2,'LineWidth',3);
set(h3,'LineWidth',3);
legend([h1 h2 h3],{'Least Squares','Huber','Student T'});
