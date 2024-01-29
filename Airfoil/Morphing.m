function [M,TFD1]=Morphing(varargin)
% digits(64)
warning('off','MATLAB:polyshape:repairedBySimplify')
w_tot=0;
for i=1:nargin
   if i==1
       M=varargin{i}*readmatrix(sprintf('./BaseShapes/%d.txt', i));
   else
       M=M+varargin{i}*readmatrix(sprintf('./BaseShapes/%d.txt', i));
   end
   w_tot=w_tot+varargin{i};
end
M=M./w_tot;

%% Intersection Reconfiguration And Stiffening (3 Loops)
TFD1 = false;
Inter_flag=0;
x=double(M(:,1));
y=double(M(:,2));

while TFD1==false

[x,y] = IntersectionReconfiguration(x,y);

if Inter_flag>5
    break;
end

x=double(x);
y=double(y);

L=[x';y'];
p = Intersection(L);

% Smoothing
y=smooth(y,3);

TFD1 = isempty(p);
Inter_flag=Inter_flag+1;

end

M(:,1)=double(x);
M(:,2)=double(y);

end

%%LOCAL SUBROUTINES
function [x,y] = IntersectionReconfiguration(x,y)

x=double(x);
y=double(y);
L=[x';y'];

p = Intersection(L);

for t=1:size(p,2)

% Redundant
if size(p,2)==0
break;
end

tx=p(1,t);
ty=p(2,t);

%% Reconfiguration
% Vector Position Start 

for r=1:size(x,1)-1

if and(or(and(x(r,1)<=tx,x(r+1,1)>=tx),and(x(r,1)>=tx,x(r+1,1)<=tx)),or(and(y(r,1)<=ty,y(r+1,1)>=ty),and(y(r,1)>=ty,y(r+1,1)<=ty)))
    break;
end

end

% Vector Position End 

for q=size(x,1):-1:2

if and(or(and(x(q,1)<=tx,x(q-1,1)>=tx),and(x(q,1)>=tx,x(q-1,1)<=tx)),or(and(y(q,1)<=ty,y(q-1,1)>=ty),and(y(q,1)>=ty,y(q-1,1)<=ty)))
    break;
end

end

q=q-1;

% Reconfiguration

rcx=(fliplr(x(r:q,:)'))';
rcy=(fliplr(y(r:q,:)'))';

x(r:q,:)=rcx;
y(r:q,:)=rcy;


%% Stiffeners

r1=r-35;
r2=r+34;
q1=q+35;
q2=q-34;

% End Intersection Control

if r2 >= q2

if r2==q2
pq2=q2-1;
q2=r2;
r2=pq2;
else
pq2=q2;
q2=r2;
r2=pq2;
end
end

% Start Intersection Control

if q1>size(x,1)
q1=1;
end
if r1<1
r1=1;    
end

% General Intersection Control

for tr=r1:r2
y(tr,1)= (((y(r2,1)-y(r1,1))/((x(r2,1)-x(r1,1))))*(x(tr,1)-x(r1,1)))+y(r1,1);
end
for tq=q1:-1:q2
y(tq,1)= (((y(q1,1)-y(q2,1))/((x(q1,1)-x(q2,1))))*(x(tq,1)-x(q1,1)))+y(q1,1);   
end

end
end

function P = Intersection(L1,varargin)
%INTERX Intersection of curves
%   P = INTERX(L1,L2) returns the intersection points of two curves L1 
%   and L2. The curves L1,L2 can be either closed or open and are described
%   by two-row-matrices, where each row contains its x- and y- coordinates.
%   The intersection of groups of curves (e.g. contour lines, multiply 
%   connected regions etc) can also be computed by separating them with a
%   column of NaNs as for example
%
%         L  = [x11 x12 x13 ... NaN x21 x22 x23 ...;
%               y11 y12 y13 ... NaN y21 y22 y23 ...]
%
%   P has the same structure as L1 and L2, and its rows correspond to the
%   x- and y- coordinates of the intersection points of L1 and L2. If no
%   intersections are found, the returned P is empty.
%
%   P = INTERX(L1) returns the self-intersection points of L1. To keep
%   the code simple, the points at which the curve is tangent to itself are
%   not included. P = INTERX(L1,L1) returns all the points of the curve 
%   together with any self-intersection points.
%   
%   Example:
%       t = linspace(0,2*pi);
%       r1 = sin(4*t)+2;  x1 = r1.*cos(t); y1 = r1.*sin(t);
%       r2 = sin(8*t)+2;  x2 = r2.*cos(t); y2 = r2.*sin(t);
%       P = InterX([x1;y1],[x2;y2]);
%       plot(x1,y1,x2,y2,P(1,:),P(2,:),'ro')

%   Author : NS
%   Version: 3.0, 21 Sept. 2010

%   Two words about the algorithm: Most of the code is self-explanatory.
%   The only trick lies in the calculation of C1 and C2. To be brief, this
%   is essentially the two-dimensional analog of the condition that needs
%   to be satisfied by a function F(x) that has a zero in the interval
%   [a,b], namely
%           F(a)*F(b) <= 0
%   C1 and C2 exactly do this for each segment of curves 1 and 2
%   respectively. If this condition is satisfied simultaneously for two
%   segments then we know that they will cross at some point. 
%   Each factor of the 'C' arrays is essentially a matrix containing 
%   the numerators of the signed distances between points of one curve
%   and line segments of the other.

    %...Argument checks and assignment of L2
    error(nargchk(1,2,nargin));
    if nargin == 1
        L2 = L1;    hF = @lt;   %...Avoid the inclusion of common points
    else
        L2 = varargin{1}; hF = @le;
    end
       
    %...Preliminary stuff
    x1  = L1(1,:)';  x2 = L2(1,:);
    y1  = L1(2,:)';  y2 = L2(2,:);
    dx1 = diff(x1); dy1 = diff(y1);
    dx2 = diff(x2); dy2 = diff(y2);
    
    %...Determine 'signed distances'   
    S1 = dx1.*y1(1:end-1) - dy1.*x1(1:end-1);
    S2 = dx2.*y2(1:end-1) - dy2.*x2(1:end-1);
    
    C1 = feval(hF,D(bsxfun(@times,dx1,y2)-bsxfun(@times,dy1,x2),S1),0);
    C2 = feval(hF,D((bsxfun(@times,y1,dx2)-bsxfun(@times,x1,dy2))',S2'),0)';

    %...Obtain the segments where an intersection is expected
    [i,j] = find(C1 & C2); 
    if isempty(i)
        P = zeros(2,0);
        return; 
    end
    
    %...Transpose and prepare for output
    i=i'; dx2=dx2'; dy2=dy2'; S2 = S2';
    L = dy2(j).*dx1(i) - dy1(i).*dx2(j);
    i = i(L~=0); j=j(L~=0); L=L(L~=0);  %...Avoid divisions by 0
    
    %...Solve system of eqs to get the common points
    P = unique([dx2(j).*S1(i) - dx1(i).*S2(j), ...
                dy2(j).*S1(i) - dy1(i).*S2(j)]./[L L],'rows')';
              
    function u = D(x,y)
        u = bsxfun(@minus,x(:,1:end-1),y).*bsxfun(@minus,x(:,2:end),y);
    end
end
