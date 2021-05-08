function q = eu2om(eu)
% EU2OM Convert from Euler angles to rotation matrix
% From Marc De Graef 3Drotations github
% INPUT: [phi1 PHI phi2] in radians
% OUTPUT: 3x3 rotation matrix

thr = 1e-10;

c1 = cos(eu(1));
c3 = cos(eu(3));
c2  = cos(eu(2));

s1 = sin(eu(1));
s3 = sin(eu(3));
s2  = sin(eu(2));

q = [ c1*c3-s1*c2*s3,  s1*c3+c1*c2*s3, s2*s3; ...
     -c1*s3-s1*c2*c3, -s1*s3+c1*c2*c3, s2*c3; ...
           s1*s2    ,      -c1*s2    ,  c2   ];

for i=1:3
  for j=1:3
    if (abs(q(i,j))< thr) 
        q(i,j) = 0.0;
    end
  end
end