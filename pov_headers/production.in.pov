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
	location <4,-90,0.5*(zhi-zlo)>
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
// randcols 0 150 0.03 0.2 0.8 1 0.8 1

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
