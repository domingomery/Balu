% ds      = Bcl_nnglm(X,d,Xt,[])  Training & Testing together
% options = Bcl_nnglm(X,d,[])     Training only
% ds      = Bcl_nnglm(Xt,options) Testing only
%
% Toolbox: Balu
%    Neural Network using a Generalized Linear Model
%
%    Design data:
%       X is a matrix with features (columns)
%       d is the ideal classification for X
%       options.method = '1, 2, 3 for 'linear','logistic' or 'softmax' (default=3)
%       options.iter is the number of iterations used in the IRLS algorithm
%       (default=10).
%
%    Test data:
%       Xt is a matrix with features (columns)
%
%    Output:
%       ds is the classification on test data
%       options.net contains information about the neural network
%       (from function Bglmtrain).
%       options.dmin contains min(d).
%       options.string is a 8 character string that describes the performed
%       classification (e.g., 'nnglm,3 ' means softmax - neural network).
%
%    Example: Training & Test together:
%       load datagauss             % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)      % plot feature space
%       op.method = 2;
%       op.iter = 12;
%       ds = Bcl_nnglm(X,d,Xt,op); % logistic - neural network
%       p = Bev_performance(ds,dt) % performance on test data
%
%    Example: Training only
%       load datagauss             % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)      % plot feature space
%       op.method = 2;
%       op.iter = 12;
%       op = Bcl_nnglm(X,d,op);    % logistic - neural network - training
%
%    Example: Testing only (after training only example):
%       ds = Bcl_nnglm(Xt,op);     % logistic - neural network - testing
%       p = Bev_performance(ds,dt) % performance on test data
%
%
%    Implementation based on NetLab Toolbox (included in this file as functions).
%
%    http://www.ncrg.aston.ac.uk/netlab/index.php
%    Copyright (c) 1996-2001, Ian T. Nabney, All rights reserved.
%
%    Nabney, I.T. (2003): Netlab: Algorithms for Pattern Recognition,
%    Advances in Pattern Recognition, Springer.
%
% D.Mery, PUC-DCC, 2010
% http://dmery.ing.puc.cl

function [ds,options] = Bcl_nnglm(varargin)
[train,test,X,d,Xt,options] = Bcl_construct(varargin{:});
options.string = sprintf('nnglm,%d ',options.method);
if train
    switch options.method
        case 1
            smethod = 'linear';
        case 2
            smethod = 'logistic';
        case 3
            smethod = 'softmax';
        otherwise
            error('Bcl_nnglm: method must be 1, 2 or 3 and not %d',method');
    end

    dmin = min(d);
    dmax = max(d);
    d1 = d-dmin+1;
    k = dmax-dmin+1;
    m = size(X,2);
    id = eye(k);
    targets = id(d1,:);
    net = Bglm(m, k, smethod);
    ops = foptions;
    ops(1) = 0;
    ops(14) = options.iter;
    options.net = Bglmtrain(net, ops, X, targets);
    options.dmin = dmin;
    ds = options;    
end
if test
    Z = Bglmfwd(options.net, Xt);
    [foo , class] = max(Z,[],2);
    ds = class+options.dmin-1;
end

end




function net = Bglm(nin, nout, outfunc, prior, beta)
%BGLM	Create a generalized linear model.
%
%	Description
%
%	NET = GLM(NIN, NOUT, FUNC) takes the number of inputs and outputs for
%	a generalized linear model, together with a string FUNC which
%	specifies the output unit activation function, and returns a data
%	structure NET. The weights are drawn from a zero mean, isotropic
%	Gaussian, with variance scaled by the fan-in of the output units.
%	This makes use of the Matlab function RANDN and so the seed for the
%	random weight initialization can be  set using RANDN('STATE', S)
%	where S is the seed value. The optional argument ALPHA sets the
%	inverse variance for the weight initialization.
%
%	The fields in NET are
%	  type = 'glm'
%	  nin = number of inputs
%	  nout = number of outputs
%	  nwts = total number of weights and biases
%	  actfn = string describing the output unit activation function:
%	      'linear'
%	      'logistic'
%	      'softmax'
%	  w1 = first-layer weight matrix
%	  b1 = first-layer bias vector
%
%	NET = Bglm(NIN, NOUT, FUNC, PRIOR), in which PRIOR is a scalar, allows
%	the field  NET.ALPHA in the data structure NET to be set,
%	corresponding  to a zero-mean isotropic Gaussian prior with inverse
%	variance with value PRIOR. Alternatively, PRIOR can consist of a data
%	structure with fields ALPHA and INDEX, allowing individual Gaussian
%	priors to be set over groups of weights in the network. Here ALPHA is
%	a column vector in which each element corresponds to a  separate
%	group of weights, which need not be mutually exclusive.  The
%	membership of the groups is defined by the matrix INDEX in which the
%	columns correspond to the elements of ALPHA. Each column has one
%	element for each weight in the matrix, in the order defined by the
%	function BglmPAK, and each element is 1 or 0 according to whether the
%	weight is a member of the corresponding group or not.
%
%	NET = Bglm(NIN, NOUT, FUNC, PRIOR, BETA) also sets the  additional
%	field NET.BETA in the data structure NET, where beta corresponds to
%	the inverse noise variance.
%
%	See also
%	BglmPAK, BglmUNPAK, BglmFWD, BglmERR, BglmGRAD, BglmTRAIN
%

%	Copyright (c) Ian T Nabney (1996-2001)

net.type = 'glm';
net.nin = nin;
net.nout = nout;
net.nwts = (nin + 1)*nout;

outtfns = {'linear', 'logistic', 'softmax'};

if sum(strcmp(outfunc, outtfns)) == 0
  error('Undefined activation function. Exiting.');
else
  net.outfn = outfunc;
end

if nargin > 3
  if isstruct(prior)
    net.alpha = prior.alpha;
    net.index = prior.index;
  elseif size(prior) == [1 1]
    net.alpha = prior;
  else
    error('prior must be a scalar or structure');
  end
end
  
net.w1 = randn(nin, nout)/sqrt(nin + 1);
net.b1 = randn(1, nout)/sqrt(nin + 1);

if nargin == 5
  net.beta = beta;
end
end

function [e, edata, eprior, y, a] = Bglmerr(net, x, t)
%BglmERR	Evaluate error function for generalized linear model.
%
%	Description
%	 E = BglmERR(NET, X, T) takes a generalized linear model data
%	structure NET together with a matrix X of input vectors and a matrix
%	T of target vectors, and evaluates the error function E. The choice
%	of error function corresponds to the output unit activation function.
%	Each row of X corresponds to one input vector and each row of T
%	corresponds to one target vector.
%
%	[E, EDATA, EPRIOR, Y, A] = BglmERR(NET, X, T) also returns the data
%	and prior components of the total error.
%
%	[E, EDATA, EPRIOR, Y, A] = BglmERR(NET, X) also returns a matrix Y
%	giving the outputs of the models and a matrix A  giving the summed
%	inputs to each output unit, where each row corresponds to one
%	pattern.
%
%	See also
%	Bglm, BglmPAK, BglmUNPAK, BglmFWD, BglmGRAD, BglmTRAIN
%

%	Copyright (c) Ian T Nabney (1996-2001)

% Check arguments for consistency
errstring = Bconsist(net, 'glm', x, t);
if ~isempty(errstring);
  error(errstring);
end

[y, a] = Bglmfwd(net, x);

switch net.outfn

  case 'linear'  	% Linear outputs
    edata = 0.5*sum(sum((y - t).^2));

  case 'logistic'  	% Logistic outputs
    edata = - sum(sum(t.*log(y) + (1 - t).*log(1 - y)));

  case 'softmax'   	% Softmax outputs
    edata = - sum(sum(t.*log(y)));

  otherwise
    error(['Unknown activation function ', net.outfn]);
end

[e, edata, eprior] = Berrbayes(net, edata);
end

function [y, a] = Bglmfwd(net, x)
%BglmFWD	Forward propagation through generalized linear model.
%
%	Description
%	Y = BglmFWD(NET, X) takes a generalized linear model data structure
%	NET together with a matrix X of input vectors, and forward propagates
%	the inputs through the network to generate a matrix Y of output
%	vectors. Each row of X corresponds to one input vector and each row
%	of Y corresponds to one output vector.
%
%	[Y, A] = BglmFWD(NET, X) also returns a matrix A  giving the summed
%	inputs to each output unit, where each row corresponds to one
%	pattern.
%
%	See also
%	Bglm, BglmPAK, BglmUNPAK, BglmERR, BglmGRAD
%

%	Copyright (c) Ian T Nabney (1996-2001)

% Check arguments for consistency
errstring = Bconsist(net, 'glm', x);
if ~isempty(errstring);
  error(errstring);
end

ndata = size(x, 1);

a = x*net.w1 + ones(ndata, 1)*net.b1;

switch net.outfn

  case 'linear'     % Linear outputs
    y = a;

  case 'logistic'   % Logistic outputs
    % Prevent overflow and underflow: use same bounds as Bglmerr
    % Ensure that log(1-y) is computable: need exp(a) > eps
    maxcut = -log(eps);
    % Ensure that log(y) is computable
    mincut = -log(1/realmin - 1);
    a = min(a, maxcut);
    a = max(a, mincut);
    y = 1./(1 + exp(-a));

  case 'softmax'   	% Softmax outputs
    nout = size(a,2);
    % Prevent overflow and underflow: use same bounds as Bglmerr
    % Ensure that sum(exp(a), 2) does not overflow
    maxcut = log(realmax) - log(nout);
    % Ensure that exp(a) > 0
    mincut = log(realmin);
    a = min(a, maxcut);
    a = max(a, mincut);
    temp = exp(a);
    y = temp./(sum(temp, 2)*ones(1,nout));
    % Ensure that log(y) is computable
    y(y<realmin) = realmin;

  otherwise
    error(['Unknown activation function ', net.outfn]);
end
end


function [g, gdata, gprior] = Bglmgrad(net, x, t)
%BglmGRAD Evaluate gradient of error function for generalized linear model.
%
%	Description
%	G = BglmGRAD(NET, X, T) takes a generalized linear model data
%	structure NET  together with a matrix X of input vectors and a matrix
%	T of target vectors, and evaluates the gradient G of the error
%	function with respect to the network weights. The error function
%	corresponds to the choice of output unit activation function. Each
%	row of X corresponds to one input vector and each row of T
%	corresponds to one target vector.
%
%	[G, GDATA, GPRIOR] = BglmGRAD(NET, X, T) also returns separately  the
%	data and prior contributions to the gradient.
%
%	See also
%	Bglm, BglmPAK, BglmUNPAK, BglmFWD, BglmERR, BglmTRAIN
%

%	Copyright (c) Ian T Nabney (1996-2001)

% Check arguments for consistency
errstring = Bconsist(net, 'glm', x, t);
if ~isempty(errstring);
  error(errstring);
end

y = Bglmfwd(net, x);
delout = y - t;

gw1 = x'*delout;
gb1 = sum(delout, 1);

gdata = [gw1(:)', gb1];

[g, gdata, gprior] = Bgbayes(net, gdata);
end

function [h, hdata] = Bglmhess(net, x, t, hdata)
%BglmHESS Evaluate the Hessian matrix for a generalised linear model.
%
%	Description
%	H = BglmHESS(NET, X, T) takes a Bglm network data structure NET,   a
%	matrix X of input values, and a matrix T of target values and returns
%	the full Hessian matrix H corresponding to the second derivatives of
%	the negative log posterior distribution, evaluated for the current
%	weight and bias values as defined by NET. Note that the target data
%	is not required in the calculation, but is included to make the
%	interface uniform with NETHESS.  For linear and logistic outputs, the
%	computation is very simple and is  done (in effect) in one line in
%	BglmTRAIN.
%
%	[H, HDATA] = BglmHESS(NET, X, T) returns both the Hessian matrix H and
%	the contribution HDATA arising from the data dependent term in the
%	Hessian.
%
%	H = BglmHESS(NET, X, T, HDATA) takes a network data structure NET, a
%	matrix X of input values, and a matrix T of  target values, together
%	with the contribution HDATA arising from the data dependent term in
%	the Hessian, and returns the full Hessian matrix H corresponding to
%	the second derivatives of the negative log posterior distribution.
%	This version saves computation time if HDATA has already been
%	evaluated for the current weight and bias values.
%
%	See also
%	Bglm, BglmTRAIN, HESSCHEK, NETHESS
%

%	Copyright (c) Ian T Nabney (1996-2001)

% Check arguments for consistency
errstring = Bconsist(net, 'glm', x, t);
if ~isempty(errstring);
  error(errstring);
end

ndata = size(x, 1);
nparams = net.nwts;
nout = net.nout;
p = Bglmfwd(net, x);
inputs = [x ones(ndata, 1)];

if nargin == 3
   hdata = zeros(nparams);	% Full Hessian matrix
   % Calculate data component of Hessian
   switch net.outfn

   case 'linear'
      % No weighting function here
      out_hess = [x ones(ndata, 1)]'*[x ones(ndata, 1)];
      for j = 1:nout
         hdata = rearrange_hess(net, j, out_hess, hdata);
      end
   case 'logistic'
      % Each output is independent
      e = ones(1, net.nin+1);
      link_deriv = p.*(1-p);
      out_hess = zeros(net.nin+1);
      for j = 1:nout
         inputs = [x ones(ndata, 1)].*(sqrt(link_deriv(:,j))*e);
         out_hess = inputs'*inputs;   % Hessian for this output
         hdata = rearrange_hess(net, j, out_hess, hdata);
      end
      
   case 'softmax'
      bb_start = nparams - nout + 1;	% Start of bias weights block
      ex_hess = zeros(nparams);	% Contribution to Hessian from single example
      for m = 1:ndata
         X = x(m,:)'*x(m,:);
         a = diag(p(m,:))-((p(m,:)')*p(m,:));
         ex_hess(1:nparams-nout,1:nparams-nout) = kron(a, X);
         ex_hess(bb_start:nparams, bb_start:nparams) = a.*ones(net.nout, net.nout);
         temp = kron(a, x(m,:));
         ex_hess(bb_start:nparams, 1:nparams-nout) = temp;
         ex_hess(1:nparams-nout, bb_start:nparams) = temp';
         hdata = hdata + ex_hess;
      end
    otherwise
      error(['Unknown activation function ', net.outfn]);
    end
end

[h, hdata] = Bhbayes(net, hdata);
end

function hdata = rearrange_hess(net, j, out_hess, hdata)

% Because all the biases come after all the input weights,
% we have to rearrange the blocks that make up the network Hessian.
% This function assumes that we are on the jth output and that all outputs
% are independent.

bb_start = net.nwts - net.nout + 1;	% Start of bias weights block
ob_start = 1+(j-1)*net.nin; 	% Start of weight block for jth output
ob_end = j*net.nin;         	% End of weight block for jth output
b_index = bb_start+(j-1);   	% Index of bias weight
% Put input weight block in right place
hdata(ob_start:ob_end, ob_start:ob_end) = out_hess(1:net.nin, 1:net.nin);
% Put second derivative of bias weight in right place
hdata(b_index, b_index) = out_hess(net.nin+1, net.nin+1);
% Put cross terms (input weight v bias weight) in right place
hdata(b_index, ob_start:ob_end) = out_hess(net.nin+1,1:net.nin);
hdata(ob_start:ob_end, b_index) = out_hess(1:net.nin, net.nin+1);

%return 
end


function w = Bglmpak(net)
%BglmPAK	Combines weights and biases into one weights vector.
%
%	Description
%	W = BglmPAK(NET) takes a network data structure NET and  combines them
%	into a single row vector W.
%
%	See also
%	Bglm, BglmUNPAK, BglmFWD, BglmERR, BglmGRAD
%

%	Copyright (c) Ian T Nabney (1996-2001)

errstring = Bconsist(net, 'glm');
if ~errstring
  error(errstring);
end

w = [net.w1(:)', net.b1];

end



function [net, options] = Bglmtrain(net, options, x, t)
%BglmTRAIN Specialised training of generalized linear model
%
%	Description
%	NET = BglmTRAIN(NET, OPTIONS, X, T) uses the iterative reweighted
%	least squares (IRLS) algorithm to set the weights in the generalized
%	linear model structure NET.  This is a more efficient alternative to
%	using BglmERR and BglmGRAD and a non-linear optimisation routine
%	through NETOPT. Note that for linear outputs, a single pass through
%	the  algorithm is all that is required, since the error function is
%	quadratic in the weights.  The algorithm also handles scalar ALPHA
%	and BETA terms.  If you want to use more complicated priors, you
%	should use general-purpose non-linear optimisation algorithms.
%
%	For logistic and softmax outputs, general priors can be handled,
%	although this requires the pseudo-inverse of the Hessian, giving up
%	the better conditioning and some of the speed advantage of the normal
%	form equations.
%
%	The error function value at the final set of weights is returned in
%	OPTIONS(8). Each row of X corresponds to one input vector and each
%	row of T corresponds to one target vector.
%
%	The optional parameters have the following interpretations.
%
%	OPTIONS(1) is set to 1 to display error values during training. If
%	OPTIONS(1) is set to 0, then only warning messages are displayed.  If
%	OPTIONS(1) is -1, then nothing is displayed.
%
%	OPTIONS(2) is a measure of the precision required for the value of
%	the weights W at the solution.
%
%	OPTIONS(3) is a measure of the precision required of the objective
%	function at the solution.  Both this and the previous condition must
%	be satisfied for termination.
%
%	OPTIONS(5) is set to 1 if an approximation to the Hessian (which
%	assumes that all outputs are independent) is used for softmax
%	outputs. With the default value of 0 the exact Hessian (which is more
%	expensive to compute) is used.
%
%	OPTIONS(14) is the maximum number of iterations for the IRLS
%	algorithm;  default 100.
%
%	See also
%	Bglm, BglmERR, BglmGRAD
%

%	Copyright (c) Ian T Nabney (1996-2001)

% Check arguments for consistency
errstring = Bconsist(net, 'glm', x, t);
if ~errstring
  error(errstring);
end

if(~options(14))
  options(14) = 100;
end

display = options(1);
% Do we need to test for termination?
test = (options(2) | options(3));

ndata = size(x, 1);
% Add a column of ones for the bias 
inputs = [x ones(ndata, 1)];

% Linear outputs are a special case as they can be found in one step
if strcmp(net.outfn, 'linear')
  if ~isfield(net, 'alpha')
    % Solve for the weights and biases using left matrix divide
    temp = inputs\t;
  elseif size(net.alpha == [1 1])
    if isfield(net, 'beta')
      beta = net.beta;
    else
      beta = 1.0;
    end
    % Use normal form equation
    hessian = beta*(inputs'*inputs) + net.alpha*eye(net.nin+1);
    temp = pinv(hessian)*(beta*(inputs'*t));  
  else
    error('Only scalar alpha allowed');
  end
  net.w1 = temp(1:net.nin, :);
  net.b1 = temp(net.nin+1, :);
  % Store error value in options vector
  options(8) = Bglmerr(net, x, t);
  return;
end

% Otherwise need to use iterative reweighted least squares
e = ones(1, net.nin+1);
for n = 1:options(14)

  switch net.outfn
    case 'logistic'
      if n == 1
        % Initialise model
        p = (t+0.5)/2;
	act = log(p./(1-p));
        wold = Bglmpak(net);
      end
      link_deriv = p.*(1-p);
      weights = sqrt(link_deriv); % sqrt of weights
      if (min(min(weights)) < eps)
        warning('ill-conditioned weights in Bglmtrain')
        return
      end
      z = act + (t-p)./link_deriv;
      if ~isfield(net, 'alpha')
         % Treat each output independently with relevant set of weights
         for j = 1:net.nout
	    indep = inputs.*(weights(:,j)*e);
	    dep = z(:,j).*weights(:,j);
	    temp = indep\dep;
	    net.w1(:,j) = temp(1:net.nin);
	    net.b1(j) = temp(net.nin+1);
         end
      else
	 gradient = Bglmgrad(net, x, t);
         Hessian = Bglmhess(net, x, t);
         deltaw = -gradient*pinv(Hessian);
         w = wold + deltaw;
         net = Bglmunpak(net, w);
      end
      [err, edata, eprior, p, act] = Bglmerr(net, x, t);
      if n == 1
        errold = err;
        wold = Bnetpak(net);
      else
        w = Bnetpak(net);
      end
    case 'softmax'
      if n == 1
        % Initialise model: ensure that row sum of p is one no matter
	% how many classes there are
        p = (t + (1/size(t, 2)))/2;
	act = log(p./(1-p));
      end
      if options(5) == 1 | n == 1
        link_deriv = p.*(1-p);
        weights = sqrt(link_deriv); % sqrt of weights
        if (min(min(weights)) < eps)
          warning('ill-conditioned weights in Bglmtrain')
          return
        end
        z = act + (t-p)./link_deriv;
        % Treat each output independently with relevant set of weights
        for j = 1:net.nout
          indep = inputs.*(weights(:,j)*e);
	  dep = z(:,j).*weights(:,j);
	  temp = indep\dep;
	  net.w1(:,j) = temp(1:net.nin);
	  net.b1(j) = temp(net.nin+1);
        end
        [err, edata, eprior, p, act] = Bglmerr(net, x, t);
        if n == 1
          errold = err;
          wold = Bnetpak(net);
        else
          w = Bnetpak(net);
        end
      else
	% Exact method of calculation after w first initialised
	% Start by working out Hessian
	Hessian = Bglmhess(net, x, t);
	gradient = Bglmgrad(net, x, t);
	% Now compute modification to weights
	deltaw = -gradient*pinv(Hessian);
	w = wold + deltaw;
	net = Bglmunpak(net, w);
	[err, edata, eprior, p] = Bglmerr(net, x, t);
    end

    otherwise
      error(['Unknown activation function ', net.outfn]);
   end
   if options(1)
     fprintf(1, 'Cycle %4d Error %11.6f\n', n, err)
   end
   % Test for termination
   % Terminate if error increases
   if err >  errold
     errold = err;
     w = wold;
     options(8) = err;
     fprintf(1, 'Error has increased: terminating\n')
     return;
   end
   if test & n > 1
     if (max(abs(w - wold)) < options(2) & abs(err-errold) < options(3))
       options(8) = err;
       return;
     else
       errold = err;
       wold = w;
     end
   end
end

options(8) = err;
if (options(1) >= 0)
  disp(Bmaxitmess);
end
end


function net = Bglmunpak(net, w)
%BglmUNPAK Separates weights vector into weight and bias matrices. 
%
%	Description
%	NET = BglmUNPAK(NET, W) takes a Bglm network data structure NET and  a
%	weight vector W, and returns a network data structure identical to
%	the input network, except that the first-layer weight matrix W1 and
%	the first-layer bias vector B1 have been set to the corresponding
%	elements of W.
%
%	See also
%	Bglm, BglmPAK, BglmFWD, BglmERR, BglmGRAD
%

%	Copyright (c) Ian T Nabney (1996-2001)

% Check arguments for consistency
errstring = Bconsist(net, 'glm');
if ~errstring
  error(errstring);
end

if net.nwts ~= length(w)
  error('Invalid weight vector length')
end

nin = net.nin;
nout = net.nout;
net.w1 = reshape(w(1:nin*nout), nin, nout);
net.b1 = reshape(w(nin*nout + 1: (nin + 1)*nout), 1, nout);
end



function w = Bnetpak(net)
%NETPAK	Combines weights and biases into one weights vector.
%
%	Description
%	W = NETPAK(NET) takes a network data structure NET and combines the
%	component weight matrices  into a single row vector W. The facility
%	to switch between these two representations for the network
%	parameters is useful, for example, in training a network by error
%	function minimization, since a single vector of parameters can be
%	handled by general-purpose optimization routines.  This function also
%	takes into account a MASK defined as a field in NET by removing any
%	weights that correspond to entries of 0 in the mask.
%
%	See also
%	NET, NETUNPAK, NETFWD, NETERR, NETGRAD
%

%	Copyright (c) Ian T Nabney (1996-2001)

pakstr = ['B' net.type, 'pak'];
w = feval(pakstr, net);
% Return masked subset of weights
if (isfield(net, 'mask'))
   w = w(logical(net.mask));
end
end


function [g, gdata, gprior] = Bgbayes(net, gdata)
%GBAYES	Evaluate gradient of Bayesian error function for network.
%
%	Description
%	G = GBAYES(NET, GDATA) takes a network data structure NET together
%	the data contribution to the error gradient for a set of inputs and
%	targets. It returns the regularised error gradient using any zero
%	mean Gaussian priors on the weights defined in NET.  In addition, if
%	a MASK is defined in NET, then the entries in G that correspond to
%	weights with a 0 in the mask are removed.
%
%	[G, GDATA, GPRIOR] = GBAYES(NET, GDATA) additionally returns the data
%	and prior components of the error.
%
%	See also
%	ERRBAYES, GLMGRAD, MLPGRAD, RBFGRAD
%

%	Copyright (c) Ian T Nabney (1996-2001)

% Evaluate the data contribution to the gradient.
if (isfield(net, 'mask'))
   gdata = gdata(logical(net.mask));
end
if isfield(net, 'beta')
  g1 = gdata*net.beta;
else
  g1 = gdata;
end

% Evaluate the prior contribution to the gradient.
if isfield(net, 'alpha')
   w = netpak(net);
   if size(net.alpha) == [1 1]
      gprior = w;
      g2 = net.alpha*gprior;
   else
      if (isfield(net, 'mask'))
         nindx_cols = size(net.index, 2);
         nmask_rows = size(find(net.mask), 1);
         index = reshape(net.index(logical(repmat(net.mask, ...
            1, nindx_cols))), nmask_rows, nindx_cols);
      else
         index = net.index;
      end
      
      ngroups = size(net.alpha, 1);
      gprior = index'.*(ones(ngroups, 1)*w);
      g2 = net.alpha'*gprior;
   end
else
  gprior = 0;
  g2 = 0;
end

g = g1 + g2;
end

function [h, hdata] = Bhbayes(net, hdata) 
%HBAYES	Evaluate Hessian of Bayesian error function for network.
%
%	Description
%	H = HBAYES(NET, HDATA) takes a network data structure NET together
%	the data contribution to the Hessian for a set of inputs and targets.
%	It returns the regularised Hessian using any zero mean Gaussian
%	priors on the weights defined in NET.  In addition, if a MASK is
%	defined in NET, then the entries in H that correspond to weights with
%	a 0 in the mask are removed.
%
%	[H, HDATA] = HBAYES(NET, HDATA) additionally returns the data
%	component of the Hessian.
%
%	See also
%	GBAYES, GLMHESS, MLPHESS, RBFHESS
%

%	Copyright (c) Ian T Nabney (1996-2001)

if (isfield(net, 'mask'))
  % Extract relevant entries in Hessian
  nmask_rows = size(find(net.mask), 1);
  hdata = reshape(hdata(logical(net.mask*(net.mask'))), ...
     nmask_rows, nmask_rows);
  nwts = nmask_rows;
else
  nwts = net.nwts;
end
if isfield(net, 'beta')
  h = net.beta*hdata;
else
  h = hdata;
end

if isfield(net, 'alpha')
  if size(net.alpha) == [1 1]
    h = h + net.alpha*eye(nwts);
  else
    if isfield(net, 'mask')
      nindx_cols = size(net.index, 2);
      index = reshape(net.index(logical(repmat(net.mask, ...
         1, nindx_cols))), nmask_rows, nindx_cols);
    else
      index = net.index;
    end
    h = h + diag(index*net.alpha);
  end 
end

end


function [e, edata, eprior] = Berrbayes(net, edata)
%ERRBAYES Evaluate Bayesian error function for network.
%
%	Description
%	E = BERRBAYES(NET, EDATA) takes a network data structure  NET together
%	the data contribution to the error for a set of inputs and targets.
%	It returns the regularised error using any zero mean Gaussian priors
%	on the weights defined in NET.
%
%	[E, EDATA, EPRIOR] = BERRBAYES(NET, X, T) additionally returns the
%	data and prior components of the error.
%
%	See also
%	GLMERR, MLPERR, RBFERR
%

%	Copyright (c) Ian T Nabney (1996-2001)

% Evaluate the data contribution to the error.
if isfield(net, 'beta')
  e1 = net.beta*edata;
else
  e1 = edata;
end

% Evaluate the prior contribution to the error.
if isfield(net, 'alpha')
   w = netpak(net);
   if size(net.alpha) == [1 1]
      eprior = 0.5*(w*w');
      e2 = eprior*net.alpha;
   else
      if (isfield(net, 'mask'))
         nindx_cols = size(net.index, 2);
         nmask_rows = size(find(net.mask), 1);
         index = reshape(net.index(logical(repmat(net.mask, ...
            1, nindx_cols))), nmask_rows, nindx_cols);
      else
         index = net.index;
      end
      eprior = 0.5*(w.^2)*index;
      e2 = eprior*net.alpha;
   end
else
  eprior = 0;
  e2 = 0;
end

e = e1 + e2;
end



function errstring = Bconsist(model, type, inputs, outputs)
%BCONSIST Check that arguments are consistent.
%
%	Description
%
%	ERRSTRING = CONSIST(NET, TYPE, INPUTS) takes a network data structure
%	NET together with a string TYPE containing the correct network type,
%	a matrix INPUTS of input vectors and checks that the data structure
%	is consistent with the other arguments.  An empty string is returned
%	if there is no error, otherwise the string contains the relevant
%	error message.  If the TYPE string is empty, then any type of network
%	is allowed.
%
%	ERRSTRING = CONSIST(NET, TYPE) takes a network data structure NET
%	together with a string TYPE containing the correct  network type, and
%	checks that the two types match.
%
%	ERRSTRING = CONSIST(NET, TYPE, INPUTS, OUTPUTS) also checks that the
%	network has the correct number of outputs, and that the number of
%	patterns in the INPUTS and OUTPUTS is the same.  The fields in NET
%	that are used are
%	  type
%	  nin
%	  nout
%
%	See also
%	MLPFWD
%

%	Copyright (c) Ian T Nabney (1996-2001)

% Assume that all is OK as default
errstring = '';

% If type string is not empty
if ~isempty(type)
  % First check that model has type field
  if ~isfield(model, 'type')
    errstring = 'Data structure does not contain type field';
    return
  end
  % Check that model has the correct type
  s = model.type;
  if ~strcmp(s, type)
    errstring = ['Model type ''', s, ''' does not match expected type ''',...
	type, ''''];
    return
  end
end

% If inputs are present, check that they have correct dimension
if nargin > 2
  if ~isfield(model, 'nin')
    errstring = 'Data structure does not contain nin field';
    return
  end

  data_nin = size(inputs, 2);
  if model.nin ~= data_nin
    errstring = ['Dimension of inputs ', num2str(data_nin), ...
	' does not match number of model inputs ', num2str(model.nin)];
    return
  end
end

% If outputs are present, check that they have correct dimension
if nargin > 3
  if ~isfield(model, 'nout')
    errstring = 'Data structure does not conatin nout field';
    return
  end
  data_nout = size(outputs, 2);
  if model.nout ~= data_nout
    errstring = ['Dimension of outputs ', num2str(data_nout), ...
	' does not match number of model outputs ', num2str(model.nout)];
    return
  end

% Also check that number of data points in inputs and outputs is the same
  num_in = size(inputs, 1);
  num_out = size(outputs, 1);
  if num_in ~= num_out
    errstring = ['Number of input patterns ', num2str(num_in), ...
	' does not match number of output patterns ', num2str(num_out)];
    return
  end
end
end




function s = Bmaxitmess()
%MAXITMESS Create a standard error message when training reaches max. iterations.
%
%	Description
%	S = MAXITMESS returns a standard string that it used by training
%	algorithms when the maximum number of iterations (as specified in
%	OPTIONS(14) is reached.
%
%	See also
%	CONJGRAD, GLMTRAIN, GMMEM, GRADDESC, GTMEM, KMEANS, OLGD, QUASINEW, SCG
%

%	Copyright (c) Ian T Nabney (1996-2001)

s = 'Maximum number of iterations has been exceeded';

end
