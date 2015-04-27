/*  detect-lines, extract lines and their width from images.
    Copyright (C) 1996-1998 Carsten Steger

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2, or (at your option)
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. */

/* 	Changes Made by R. Balasubramanian for incorporating the the detect lines code to incorporate
   	within GRASP (May 10th 1999) */

/*	Translated into an ImageJ java plugin by Thorsten Wagner (Dez. 2014) */

package de.biomedical_imaging.ij.steger;

public class Threshold {
	
	static void threshold(byte[] image, long min,long width,long height,Region out)

	{
	  long   grey;
	  long   r,c,l,num,num_max;
	  boolean   inside;
	  Chord[] rl;

	  inside = false;
	  num = 0;
	  num_max = LinesUtil.INITIAL_SIZE;
	  rl =new Chord[(int) num_max]; 
	  for(int i = 0; i < num_max; i++){
		  rl[i] = new Chord();
	  }
	  out.rl = null;
	  out.num = 0;

	  for (r=0; r<height; r++) {
	    for (c=0; c<width; c++) {
	      l = LinesUtil.LINCOOR(r,c,width);
	      grey = image[(int) l];
	      if (grey >= min) {
	        if (!inside) {
	          inside = true;
	          rl[(int) num].r = (short) r;
	          rl[(int) num].cb = (short) c;
	        }
	      } else {
	        if (inside) {
	          inside = false;
	          rl[(int) num].ce = (short) (c - 1);
	          num++;
	          if (num >= num_max) {
	            num_max = (long) Math.floor((double)(num_max*LinesUtil.REALLOC_FACTOR));
	            Chord[] rlh = new Chord[(int) num_max];
	            for(int i = 0; i < rlh.length; i++){
	            	if(i < rl.length)
		        		rlh[i] = rl[i];
		        	else
		        		rlh[i] = new Chord();
	      	  	}
	            rl=rlh;
	          }
	        }
	      }
	    }
	    if (inside) {
	      inside = false;
	      rl[(int) num].ce = (short) (width-1);
	      num++;
	      if (num >= num_max) {
	        num_max = (long) Math.floor((double)(num_max*LinesUtil.REALLOC_FACTOR));
	        Chord[] rlh = new Chord[(int) num_max];
	        for(int i = 0; i < rlh.length; i++){
	        	if(i < rl.length)
	        		rlh[i] = rl[i];
	        	else
	        		rlh[i] = new Chord();
	        }
	        rl = rlh;
	      }
	    }
	  }
	  out.rl = new Chord[(int) num];
	  for(int i = 0; i < (int)num; i++){
		  out.rl[i] = rl[i];
	  }
	  out.num = num;
	}

}
