--------------------------------
-- 2.9; take f as below
-- generate I = (f,fx,fy);
-- find complex dimension CC[x,y]/I and 
-- if x5 == y5 mod I
--------------------------------
R = CC[x,y,z];
f = x^4 + x^2*y^2+y^3-x^3
fx = diff(x,f)
fy = diff(y,f)
I = ideal(f,fx,fy)
J = gens gb I
-- | y2 x2 | -- hence CC[x,y]/I has basis (1,x,y,xy) so 4-dimensional;
-- x5 == y5 mod I => true

--------------------------------------------------
-- Ch 2.10 Groebner for <x_3-x_1^5,x_2-x_1^3>
-- computational difference using Lex & GRevLex
-- normal form of x_1x_2x_3 can differ
--------------------------------------------------
RLex = QQ[x_1, x_2, x_3, MonomialOrder => Lex]
describe RLex
ILex = ideal(x_3-x_1^5,x_2-x_1^3)
gbTrace = 4
JLex = gens gb ILex
-- | x_2^5-x_3^3 x_1x_3-x_2^2 x_1x_2^3-x_3^2 x_1^2x_2-x_3 x_1^3-x_2 |
-- 8 spairs done
RGRevLex = QQ[y_1, y_2, y_3, MonomialOrder => GRevLex]
describe RGRevLex
IGRevLex = ideal(y_3-y_1^5,y_2-y_1^3)
JRevLex = gens gb IGRevLex
-- | y_2^2-y_1y_3 y_1^2y_2-y_3 y_1^3-y_2 |
-- spairs done = 4
x_1*x_2*x_3 % ILex 
-- > the normal form is x_2^3
y_1*y_2*y_3 % IGRevLex
--> the normal form is the monomial itself