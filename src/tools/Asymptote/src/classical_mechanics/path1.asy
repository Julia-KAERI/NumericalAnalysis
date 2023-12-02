import graph;
import geometry;
import math;
import settings;
import fontsize;


settings.outformat = "pdf";
defaultpen(fontsize(17pt));
defaultpen(1);

size(400,300,IgnoreAspect);

real a=1, b=1.5, c=2;
real ya=1.0, yb=2.4, yc=3.0;
real yb2=2.6, yb3=2.2;
int t=2;

real f(real x) {
    real r1=(x-b)*(x-c)/(a-b)/(a-c)*ya ;
    r1 += (x-a)*(x-c)/(b-a)/(b-c)*yb;
    r1 += (x-a)*(x-b)/(c-a)/(c-b)*yc;
    return r1;
    }

real f2(real x) {
    real r1=(x-b)*(x-c)/(a-b)/(a-c)*ya ;
    r1 += (x-a)*(x-c)/(b-a)/(b-c)*yb2;
    r1 += (x-a)*(x-b)/(c-a)/(c-b)*yc;
    return r1;
    }
real f3(real x) {
    real r1=(x-b)*(x-c)/(a-b)/(a-c)*ya ;
    r1 += (x-a)*(x-c)/(b-a)/(b-c)*yb3;
    r1 += (x-a)*(x-b)/(c-a)/(c-b)*yc;
    return r1;
    }

pair F(real x) {return (x,f(x));}

dotfactor=7;

xaxis("$t$", xmin=0.5, xmax=2.5, Arrow, ticks=Ticks(DefaultFormat,
                                        new real[] {1, 2}));
yaxis("$P$", XEquals(0.7), ymin=-0.2, ymax=4, Arrow);

// for(int i=a; i < c; ++i)
//   draw(F(i)--F(i+1));

path g=graph(f,a,c);
path g2=graph(f2,a,c);
path g3=graph(f3,a,c);

draw(g,black);   
draw(g2,blue+dashed);  
draw(g3,dashed+red);  

int n=2;

// perpendicular((t+1,f(t)),NW,blue);

// draw(F(t)--(t+1,f(t)),blue);
// draw((t+1,f(t))--F(t+1),heavygreen);

// dot(F(t));
// dot(F(t+1));


// real sx=0.08;
// real sy=1.5;
// draw("$x(t+h)-x(t)$",F(t)-sy*I--(t+1,f(t)-sy),blue+dashed,Arrows,Bars,
//      PenMargins);
// draw("$y(t+h)-y(t)$",(t+1+sx,f(t))--F(t+1)+sx,heavygreen+dashed,Arrows,Bars,PenMargins);

// labelx("$1$",a);
// labelx("$2$",c);


dot(Label("$P_1$",align=W), F(a));
dot(Label("$P_2$",align=E), F(c));
// label("$s(a)=0$",F(a),W,red);
// label("$s(t)$",F(t),N,red);
// label("$s(t+h)$",F(t+1),NW,red);
// label("$s(b)=L$",F(b),ESE,red);
// label("$P=(x(t),y(t))$",F(t),2.5W);
// label("$Q=(x(t+h),y(t+h))$",F(t+1),NE);