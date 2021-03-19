// W=1200 H=900
#version 3.7;

global_settings {
	max_trace_level 64
}

#declare xlo=0;
#declare ylo=0;
#declare zlo=0;
#declare xhi=1;
#declare yhi=1;
#declare zhi=1;

camera {
	location <4,-20,0.5*(zhi-zlo)>
	sky z
	right -0.066*x*image_width/image_height
	up 0.066*z
	look_at <0.5,0.5,0.5*(zhi-zlo)>
}

// White background
background{rgb 1}

light_source{<-8,-20,30> color rgb <0.77,0.75,0.75> shadowless}
light_source{<25,-12,12> color rgb <0.38,0.40,0.40> shadowless}
#declare f0=finish{specular 0.5 ambient 0.25}

// The following line gets parsed by random_colors.pl
// to produce mesh of various colors depending on object ID
// Second integer = the number of distinct objects
#declare t_msh0=texture{pigment{rgb <1.000,0.391,0.179>} finish{f0}}
#declare t_msh1=texture{pigment{rgb <0.847,0.437,0.007>} finish{f0}}
#declare t_msh2=texture{pigment{rgb <0.845,0.235,0.083>} finish{f0}}
#declare t_msh3=texture{pigment{rgb <0.837,0.458,0.139>} finish{f0}}
#declare t_msh4=texture{pigment{rgb <0.823,0.659,0.002>} finish{f0}}
#declare t_msh5=texture{pigment{rgb <0.844,0.486,0.136>} finish{f0}}
#declare t_msh6=texture{pigment{rgb <0.948,0.739,0.183>} finish{f0}}
#declare t_msh7=texture{pigment{rgb <0.943,0.755,0.178>} finish{f0}}
#declare t_msh8=texture{pigment{rgb <0.997,0.622,0.143>} finish{f0}}
#declare t_msh9=texture{pigment{rgb <0.879,0.310,0.125>} finish{f0}}
#declare t_msh10=texture{pigment{rgb <0.900,0.332,0.163>} finish{f0}}
#declare t_msh11=texture{pigment{rgb <0.946,0.837,0.023>} finish{f0}}
#declare t_msh12=texture{pigment{rgb <0.856,0.737,0.041>} finish{f0}}
#declare t_msh13=texture{pigment{rgb <0.908,0.834,0.153>} finish{f0}}
#declare t_msh14=texture{pigment{rgb <0.861,0.856,0.124>} finish{f0}}
#declare t_msh15=texture{pigment{rgb <0.770,0.810,0.030>} finish{f0}}
#declare t_msh16=texture{pigment{rgb <0.941,0.304,0.061>} finish{f0}}
#declare t_msh17=texture{pigment{rgb <0.884,0.402,0.143>} finish{f0}}
#declare t_msh18=texture{pigment{rgb <0.910,0.570,0.057>} finish{f0}}
#declare t_msh19=texture{pigment{rgb <0.876,0.515,0.113>} finish{f0}}
#declare t_msh20=texture{pigment{rgb <0.925,0.617,0.062>} finish{f0}}
#declare t_msh21=texture{pigment{rgb <0.877,0.713,0.079>} finish{f0}}
#declare t_msh22=texture{pigment{rgb <0.890,0.564,0.137>} finish{f0}}
#declare t_msh23=texture{pigment{rgb <0.800,0.500,0.070>} finish{f0}}
#declare t_msh24=texture{pigment{rgb <0.942,0.961,0.060>} finish{f0}}
#declare t_msh25=texture{pigment{rgb <0.905,0.729,0.055>} finish{f0}}
#declare t_msh26=texture{pigment{rgb <0.931,0.794,0.132>} finish{f0}}
#declare t_msh27=texture{pigment{rgb <0.779,0.854,0.058>} finish{f0}}
#declare t_msh28=texture{pigment{rgb <0.925,0.953,0.084>} finish{f0}}
#declare t_msh29=texture{pigment{rgb <0.820,0.987,0.143>} finish{f0}}
#declare t_msh30=texture{pigment{rgb <0.942,0.313,0.132>} finish{f0}}
#declare t_msh31=texture{pigment{rgb <0.841,0.828,0.020>} finish{f0}}
#declare t_msh32=texture{pigment{rgb <0.857,0.217,0.015>} finish{f0}}
#declare t_msh33=texture{pigment{rgb <0.936,0.521,0.117>} finish{f0}}
#declare t_msh34=texture{pigment{rgb <0.967,0.326,0.155>} finish{f0}}
#declare t_msh35=texture{pigment{rgb <0.906,0.970,0.063>} finish{f0}}
#declare t_msh36=texture{pigment{rgb <0.812,0.520,0.049>} finish{f0}}
#declare t_msh37=texture{pigment{rgb <0.830,0.905,0.172>} finish{f0}}
#declare t_msh38=texture{pigment{rgb <0.795,0.925,0.018>} finish{f0}}
#declare t_msh39=texture{pigment{rgb <0.955,0.450,0.088>} finish{f0}}
#declare t_msh40=texture{pigment{rgb <0.973,0.617,0.082>} finish{f0}}
#declare t_msh41=texture{pigment{rgb <0.890,0.329,0.148>} finish{f0}}
#declare t_msh42=texture{pigment{rgb <0.985,0.887,0.012>} finish{f0}}
#declare t_msh43=texture{pigment{rgb <0.970,0.699,0.143>} finish{f0}}
#declare t_msh44=texture{pigment{rgb <0.827,0.587,0.162>} finish{f0}}
#declare t_msh45=texture{pigment{rgb <0.897,0.511,0.018>} finish{f0}}
#declare t_msh46=texture{pigment{rgb <0.837,0.451,0.138>} finish{f0}}
#declare t_msh47=texture{pigment{rgb <0.825,0.408,0.118>} finish{f0}}
#declare t_msh48=texture{pigment{rgb <0.968,0.731,0.085>} finish{f0}}
#declare t_msh49=texture{pigment{rgb <0.904,0.979,0.163>} finish{f0}}
#declare t_msh50=texture{pigment{rgb <0.860,0.439,0.048>} finish{f0}}
#declare t_msh51=texture{pigment{rgb <0.999,0.750,0.072>} finish{f0}}
#declare t_msh52=texture{pigment{rgb <0.832,0.836,0.091>} finish{f0}}
#declare t_msh53=texture{pigment{rgb <0.991,0.416,0.188>} finish{f0}}
#declare t_msh54=texture{pigment{rgb <0.948,0.674,0.162>} finish{f0}}
#declare t_msh55=texture{pigment{rgb <0.925,0.322,0.159>} finish{f0}}
#declare t_msh56=texture{pigment{rgb <0.821,0.494,0.100>} finish{f0}}
#declare t_msh57=texture{pigment{rgb <0.801,0.754,0.084>} finish{f0}}
#declare t_msh58=texture{pigment{rgb <0.995,0.870,0.001>} finish{f0}}
#declare t_msh59=texture{pigment{rgb <0.907,0.986,0.190>} finish{f0}}
#declare t_msh60=texture{pigment{rgb <0.851,0.280,0.093>} finish{f0}}
#declare t_msh61=texture{pigment{rgb <0.999,0.458,0.069>} finish{f0}}
#declare t_msh62=texture{pigment{rgb <0.793,0.860,0.055>} finish{f0}}
#declare t_msh63=texture{pigment{rgb <0.817,0.289,0.106>} finish{f0}}
#declare t_msh64=texture{pigment{rgb <0.974,0.215,0.028>} finish{f0}}
#declare t_msh65=texture{pigment{rgb <0.963,0.406,0.084>} finish{f0}}
#declare t_msh66=texture{pigment{rgb <0.816,0.412,0.048>} finish{f0}}
#declare t_msh67=texture{pigment{rgb <0.926,0.372,0.162>} finish{f0}}
#declare t_msh68=texture{pigment{rgb <0.897,0.246,0.014>} finish{f0}}
#declare t_msh69=texture{pigment{rgb <0.997,0.983,0.135>} finish{f0}}
#declare t_msh70=texture{pigment{rgb <0.915,0.503,0.056>} finish{f0}}
#declare t_msh71=texture{pigment{rgb <0.932,0.809,0.056>} finish{f0}}
#declare t_msh72=texture{pigment{rgb <0.837,0.306,0.101>} finish{f0}}
#declare t_msh73=texture{pigment{rgb <0.943,0.951,0.046>} finish{f0}}
#declare t_msh74=texture{pigment{rgb <0.817,0.569,0.084>} finish{f0}}
#declare t_msh75=texture{pigment{rgb <0.950,0.937,0.187>} finish{f0}}
#declare t_msh76=texture{pigment{rgb <0.812,0.403,0.059>} finish{f0}}
#declare t_msh77=texture{pigment{rgb <0.893,0.272,0.131>} finish{f0}}
#declare t_msh78=texture{pigment{rgb <0.958,0.492,0.046>} finish{f0}}
#declare t_msh79=texture{pigment{rgb <0.974,0.623,0.146>} finish{f0}}
#declare t_msh80=texture{pigment{rgb <0.986,0.907,0.182>} finish{f0}}
#declare t_msh81=texture{pigment{rgb <0.858,0.305,0.151>} finish{f0}}
#declare t_msh82=texture{pigment{rgb <0.973,0.929,0.173>} finish{f0}}
#declare t_msh83=texture{pigment{rgb <0.839,0.700,0.043>} finish{f0}}
#declare t_msh84=texture{pigment{rgb <0.707,0.853,0.075>} finish{f0}}
#declare t_msh85=texture{pigment{rgb <0.847,0.942,0.017>} finish{f0}}
#declare t_msh86=texture{pigment{rgb <0.817,0.515,0.032>} finish{f0}}
#declare t_msh87=texture{pigment{rgb <0.868,0.726,0.163>} finish{f0}}
#declare t_msh88=texture{pigment{rgb <0.942,0.826,0.106>} finish{f0}}
#declare t_msh89=texture{pigment{rgb <0.981,0.325,0.052>} finish{f0}}
#declare t_msh90=texture{pigment{rgb <0.896,0.648,0.076>} finish{f0}}
#declare t_msh91=texture{pigment{rgb <0.934,0.830,0.050>} finish{f0}}
#declare t_msh92=texture{pigment{rgb <0.993,0.386,0.129>} finish{f0}}
#declare t_msh93=texture{pigment{rgb <0.884,0.587,0.152>} finish{f0}}
#declare t_msh94=texture{pigment{rgb <0.929,0.765,0.108>} finish{f0}}
#declare t_msh95=texture{pigment{rgb <0.997,0.271,0.073>} finish{f0}}
#declare t_msh96=texture{pigment{rgb <0.886,0.694,0.068>} finish{f0}}
#declare t_msh97=texture{pigment{rgb <0.948,0.778,0.083>} finish{f0}}
#declare t_msh98=texture{pigment{rgb <0.931,0.319,0.159>} finish{f0}}
#declare t_msh99=texture{pigment{rgb <0.988,0.461,0.091>} finish{f0}}
#declare t_msh100=texture{pigment{rgb <0.837,0.285,0.126>} finish{f0}}
#declare t_msh101=texture{pigment{rgb <0.920,0.678,0.040>} finish{f0}}
#declare t_msh102=texture{pigment{rgb <0.929,0.384,0.017>} finish{f0}}
#declare t_msh103=texture{pigment{rgb <0.918,0.410,0.120>} finish{f0}}
#declare t_msh104=texture{pigment{rgb <0.810,0.587,0.090>} finish{f0}}
#declare t_msh105=texture{pigment{rgb <0.956,0.816,0.018>} finish{f0}}
#declare t_msh106=texture{pigment{rgb <0.869,0.576,0.145>} finish{f0}}
#declare t_msh107=texture{pigment{rgb <0.843,0.237,0.045>} finish{f0}}
#declare t_msh108=texture{pigment{rgb <0.881,0.361,0.173>} finish{f0}}
#declare t_msh109=texture{pigment{rgb <0.907,0.619,0.052>} finish{f0}}
#declare t_msh110=texture{pigment{rgb <0.819,0.605,0.106>} finish{f0}}
#declare t_msh111=texture{pigment{rgb <0.994,0.859,0.147>} finish{f0}}
#declare t_msh112=texture{pigment{rgb <0.873,0.422,0.165>} finish{f0}}
#declare t_msh113=texture{pigment{rgb <0.796,0.914,0.100>} finish{f0}}
#declare t_msh114=texture{pigment{rgb <0.816,0.354,0.073>} finish{f0}}
#declare t_msh115=texture{pigment{rgb <0.837,0.829,0.057>} finish{f0}}
#declare t_msh116=texture{pigment{rgb <0.906,0.268,0.103>} finish{f0}}
#declare t_msh117=texture{pigment{rgb <0.817,0.684,0.043>} finish{f0}}
#declare t_msh118=texture{pigment{rgb <0.996,0.504,0.069>} finish{f0}}
#declare t_msh119=texture{pigment{rgb <0.851,0.872,0.157>} finish{f0}}
#declare t_msh120=texture{pigment{rgb <0.830,0.844,0.121>} finish{f0}}
#declare t_msh121=texture{pigment{rgb <0.924,0.923,0.171>} finish{f0}}
#declare t_msh122=texture{pigment{rgb <0.948,0.388,0.099>} finish{f0}}
#declare t_msh123=texture{pigment{rgb <0.830,0.637,0.147>} finish{f0}}
#declare t_msh124=texture{pigment{rgb <0.863,0.840,0.056>} finish{f0}}
#declare t_msh125=texture{pigment{rgb <0.816,0.389,0.032>} finish{f0}}
#declare t_msh126=texture{pigment{rgb <0.886,0.684,0.001>} finish{f0}}
#declare t_msh127=texture{pigment{rgb <0.815,0.517,0.121>} finish{f0}}
#declare t_msh128=texture{pigment{rgb <0.997,0.605,0.192>} finish{f0}}
#declare t_msh129=texture{pigment{rgb <0.884,0.590,0.051>} finish{f0}}
#declare t_msh130=texture{pigment{rgb <0.849,0.849,0.165>} finish{f0}}
#declare t_msh131=texture{pigment{rgb <0.805,0.417,0.011>} finish{f0}}
#declare t_msh132=texture{pigment{rgb <0.872,0.327,0.038>} finish{f0}}
#declare t_msh133=texture{pigment{rgb <0.938,0.918,0.148>} finish{f0}}
#declare t_msh134=texture{pigment{rgb <0.809,0.900,0.015>} finish{f0}}
#declare t_msh135=texture{pigment{rgb <0.820,0.324,0.148>} finish{f0}}
#declare t_msh136=texture{pigment{rgb <0.888,0.468,0.020>} finish{f0}}
#declare t_msh137=texture{pigment{rgb <0.978,0.339,0.059>} finish{f0}}
#declare t_msh138=texture{pigment{rgb <0.917,0.968,0.087>} finish{f0}}
#declare t_msh139=texture{pigment{rgb <0.995,0.905,0.067>} finish{f0}}
#declare t_msh140=texture{pigment{rgb <0.886,0.733,0.052>} finish{f0}}
#declare t_msh141=texture{pigment{rgb <0.935,0.660,0.153>} finish{f0}}
#declare t_msh142=texture{pigment{rgb <0.998,0.455,0.034>} finish{f0}}
#declare t_msh143=texture{pigment{rgb <0.987,0.499,0.042>} finish{f0}}
#declare t_msh144=texture{pigment{rgb <0.946,0.812,0.013>} finish{f0}}
#declare t_msh145=texture{pigment{rgb <0.913,0.713,0.151>} finish{f0}}
#declare t_msh146=texture{pigment{rgb <0.687,0.814,0.045>} finish{f0}}
#declare t_msh147=texture{pigment{rgb <0.748,0.820,0.090>} finish{f0}}
#declare t_msh148=texture{pigment{rgb <0.856,0.214,0.009>} finish{f0}}
#declare t_msh149=texture{pigment{rgb <0.889,0.417,0.176>} finish{f0}}
#declare t_msh150=texture{pigment{rgb <0.918,0.961,0.013>} finish{f0}}

#include "msh.pov"

CYL:#declare s=0.0015;
CYL:
CYL:#include "cyl.pov"
CYL:
CYL:union{cyls0
CYL:	texture {
CYL:		pigment{rgb <0.2,0.3,0.8>}
CYL:		finish{specular 0.4 ambient 0.4}
CYL:	}
CYL:}

TRA:#declare r=0.005;
TRA:#declare cr=0.006;
TRA:#declare csr=0.006;
TRA:#declare tr=0.5;
TRA:union {
TRA:#include "sph.pov"
TRA:	texture {
TRA:		pigment{
TRA:			gradient y
TRA:			color_map {
TRA:				[0 color rgbt <1,1,0.7,tr>]
TRA:				[0.2 color rgbt <1,1,0.3,tr>]
TRA:				[0.55 color rgbt <0.8,0.4,0.1,tr>]
TRA:				[0.9 color rgbt <0.9,0.15,0.05,tr>]
TRA:				[1 color rgbt <0.7,0.1,0.6,tr>]
TRA:			}
TRA:			scale 1+2*r
TRA:			translate -r
TRA:		}
TRA:		finish{specular 0.4 ambient 0.3}
TRA:	}
TRA:}

#declare ccr=0.006;
#declare ccsr=0.006;
#declare boxr=0.004;
//  box{
//      <0,0,0>,
//      <1,boxr,1>
//      texture {
//          pigment{
//              rgb <0.69,0.69,0.69>
//          }
//          finish{specular 1 ambient 0.2 phong 0.7 reflection 0.1}
//      }
//  }
box{
    <xlo,yhi-boxr,zlo>,
    <xhi, yhi, zhi>
    texture {
        pigment{
            rgb <0.69,0.69,0.69>
        }
        finish{specular 1 ambient 0.2 phong 0.7 reflection 0.1}
    }
}
box{
    <xlo, ylo, zlo>,
    <xhi, yhi, boxr>
    texture {
        pigment{
            rgb <0.69,0.69,0.69>
        }
        finish{specular 1 ambient 0.2 phong 0.7 reflection 0.1}
    }
}
box{
    <xlo, ylo, zlo>,
    <boxr,yhi,zhi>
    texture {
        pigment{
            rgb <0.69,0.69,0.69>
        }
        finish{specular 1 ambient 0.2 phong 0.7 reflection 0.1}
    }
}
union {
	cylinder {<xlo,ylo,zlo>,<xhi,ylo,zlo>,ccr}
	cylinder {<xlo,ylo,zlo>,<xlo,yhi,zlo>,ccr}
	cylinder {<xlo,ylo,zlo>,<xlo,ylo,zhi>,ccr}
	cylinder {<xhi,yhi,zlo>,<xhi,yhi,zhi>,ccr}
	cylinder {<xhi,ylo,zhi>,<xhi,yhi,zhi>,ccr}
	cylinder {<xlo,yhi,zhi>,<xhi,yhi,zhi>,ccr}
	cylinder {<xhi,ylo,zlo>,<xhi,yhi,zlo>,ccr}
	cylinder {<xhi,ylo,zlo>,<xhi,ylo,zhi>,ccr}
	cylinder {<xlo,yhi,zlo>,<xhi,yhi,zlo>,ccr}
	cylinder {<xlo,yhi,zlo>,<xlo,yhi,zhi>,ccr}
	cylinder {<xlo,ylo,zhi>,<xhi,ylo,zhi>,ccr}
	cylinder {<xlo,ylo,zhi>,<xlo,yhi,zhi>,ccr}
	sphere {<xlo,ylo,zlo>,ccsr}
	sphere {<xlo,ylo,zhi>,ccsr}
	sphere {<xlo,yhi,zlo>,ccsr}
	sphere {<xlo,yhi,zhi>,ccsr}
	sphere {<xhi,ylo,zlo>,ccsr}
	sphere {<xhi,ylo,zhi>,ccsr}
	sphere {<xhi,yhi,zlo>,ccsr}
	sphere {<xhi,yhi,zhi>,ccsr}
	texture {
		pigment{
			rgb <0,0,0>
		}
		finish{specular 0.3 phong 0.7 ambient 0.2 reflection 0.1}
	}
}
