function ltid_mobt_test

z=tf('z');
tf(ltid_mobt(1+1/z,1/5))
fprintf('\n     ... should be (0.8z+0.8)/(z-0.2)\n')
tf(ltid_mobt(1/z^3,-1/2))
fprintf('\n     ... should be (0.125z^3+0.75z^2+1.5z+1)/(z^3+1.5z62+0.75z+0.125)\n')