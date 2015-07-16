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
