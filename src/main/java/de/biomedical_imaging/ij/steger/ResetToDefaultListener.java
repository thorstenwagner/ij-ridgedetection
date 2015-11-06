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
    aint with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. */

/* 	Changes Made by R. Balasubramanian for incorporating the the detect lines code to incorporate
   	within GRASP (May 10th 1999) */

/*	Translated into an ImageJ java plugin by Thorsten Wagner (Dez. 2014) */

package de.biomedical_imaging.ij.steger;

import ij.IJ;

import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

public class ResetToDefaultListener implements ActionListener {

	GenericDialogPlus gd;
	public ResetToDefaultListener(GenericDialogPlus gd) {
		this.gd = gd;
	}
	@Override
	public void actionPerformed(ActionEvent arg0) {
		
		//Set settings to default
		
		TextField textLineWidth = (TextField) gd.getNumericFields().get(0);
		textLineWidth.setText("" + IJ.d2s(Lines_.lineWidthDefault, 2));
		textLineWidth.setEditable(true);
		
		TextField textHighCon = (TextField) gd.getNumericFields().get(1);
		textHighCon.setText("" + IJ.d2s(Lines_.contrastHighDefault, 0));
		textHighCon.setEditable(true);
		
		TextField textLowCon = (TextField) gd.getNumericFields().get(2);
		textLowCon.setText("" + IJ.d2s(Lines_.contrastLowDefault, 0));
		textLowCon.setEditable(true);
		
		TextField textSigma = (TextField) gd.getNumericFields().get(3);
		textSigma.setText("" + IJ.d2s(Lines_.sigmaDefault, 2));
		textSigma.setEditable(true);
		
		TextField textLowThresh = (TextField) gd.getNumericFields().get(4);
		textLowThresh.setText("" + IJ.d2s(Lines_.lowerThreshDefault, 2));
		textLowThresh.setEditable(true);
		
		TextField textUppThresh = (TextField) gd.getNumericFields().get(5);
		textUppThresh.setText("" + IJ.d2s(Lines_.upperThreshDefault, 2));
		textUppThresh.setEditable(true);
		
	}

}
