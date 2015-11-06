/*
 * #%L
 * Ridge Detection plugin for ImageJ
 * %%
 * Copyright (C) 2014 - 2015 Thorsten Wagner (ImageJ java plugin), 1996-1998 Carsten Steger (original C code), 1999 R. Balasubramanian (detect lines code to incorporate within GRASP)
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 2 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
 * #L%
 */
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
