function [X_est, errPlot] = fista(b,f,fT,step,sizeX,opts)
%% Author:
% Vivek Boominathan
% Rice University
% vivekb@rice.edu


t_prv = 1;
X_prv = zeros([sizeX, 1],'like',b);
V_prv = X_prv;

k = 1;
fidAcc = [];
kAcc = [];

fidFn = @(X) sum(reshape(f(X) - b,[],1).^2) + opts.lam2*sum(X(:).^2);
gradFn = @(X) - 2*(fT(b)) + 2*fT(f(X)) + 2*opts.lam2*X;

%% Reconstruction
fidFn_nxt = fidFn(X_prv);

if opts.showFigs
    fh1 = figure;
    fh2 = figure;
end
errPlot = [];
tic,
while 1
    fidFn_prv = fidFn_nxt;

    X_nxt = V_prv - step * gradFn(V_prv);
    X_nxt = shrinkOp(X_nxt, step*opts.lam1);
    if (opts.nonneg)
        X_nxt = max(X_nxt, 0);
    end

    % Neterov's
    t_nxt = 0.5*(1+sqrt(1+4*t_prv^2));
    V_nxt = X_nxt + (t_prv-1)*(X_nxt-X_prv)/t_nxt;

    fidFn_nxt = fidFn(X_nxt);
    if fidFn_nxt > fidFn_prv
        t_prv = 1;
    else
        t_prv = t_nxt;
    end

    X_prv = X_nxt;
    V_prv = V_nxt;

    k = k+1;

    fidAcc = [fidAcc fidFn_nxt];
    kAcc = [kAcc k];
    frames = 5; %iterations befpre figure update
    if mod(k,frames) == 0 && opts.showFigs
        fprintf('Iter: %d \n',k);
        X_out = gather(max(X_prv,[],3));
        scaleX = max(X_out(:)); %Scaling
        X_out = X_out/scaleX; %Scaling
        figure(fh1), imshow(X_out), title('Max proj Recon'); drawnow;          
        figure(fh2), semilogy(kAcc, fidAcc), title('Error'); drawnow;
    end

    if k>opts.maxItr, break; end
    
    relErr = abs(fidFn_nxt-fidFn_prv)/fidFn_prv;
    errPlot = [errPlot fidFn_nxt];
    if relErr < opts.tol
        break;
    end 
end

X_est = X_prv;
        
end

function [Ev,Em] = eigMax_AtA(A,At,maxit,lmbd2)

szX = [size(A,2),1];

X = 10*randn(szX);

for ii = 1:maxit
    Xnew = At*(A*X) + lmbd2*X;
    Em = (Xnew(:)'*X(:))/(X(:)'*X(:));
    X = Xnew/norm(Xnew(:));
end
Ev = X;

end

function s = shrinkOp(z,lmbd,nonneg)
s = max(abs(z)-lmbd,0).*sign(z);
end