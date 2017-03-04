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

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.apache.commons.lang3.mutable.MutableInt;

import ij.IJ;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class LineDetector {
	boolean isDarkLine = false;
	boolean doCorrectPosition = false;
	boolean doEstimateWidth = false;
	boolean doExtendLine = false;
	private Options opts = null;
	private Junctions junctions;
	private Lines lines;
	Set<Integer> alreadyProcessedJunctionPoints;
	boolean bechatty = false;

	/**
	 * 
	 * @param ip
	 *            The image where the lines should be detected
	 * @param sigma
	 *            A value which depends on the line width: sigma >=
	 *            width/(2*sqrt(3))
	 * @param upperThresh
	 *            Upper hysteresis thresholds used in the linking algorithm
	 *            (Depends on the maximum contour brightness (greyvalue))
	 * @param lowerThresh
	 *            Lower hysteresis thresholds used in the linking algorithm
	 *            (Depends on the minimum contour brightness (greyvalue)).
	 * @param isDarkLine
	 *            True if the line darker than the background
	 * @param doCorrectPosition
	 *            Determines whether the line width and position correction
	 *            should be applied
	 * @param doEstimateWidth
	 *            Determines whether the line width should be extracted
	 * @param doExtendLine
	 *            Extends the detect lines to find more junction points
	 * @return An arraylist with line objects
	 */
	public Lines detectLines(ImageProcessor ip, double sigma,
			double upperThresh, double lowerThresh, double minLength, double maxLength, boolean isDarkLine,
			boolean doCorrectPosition, boolean doEstimateWidth,
			boolean doExtendLine) {
		return detectLines(ip, sigma, upperThresh, lowerThresh, minLength, maxLength,isDarkLine,
			doCorrectPosition, doEstimateWidth, doExtendLine, OverlapOption.NONE);
	}

	public Lines detectLines(ImageProcessor ip, double sigma,
		double upperThresh, double lowerThresh, double minLength, double maxLength, boolean isDarkLine,
		boolean doCorrectPosition, boolean doEstimateWidth,
		boolean doExtendLine, OverlapOption overlapOption) {
	this.isDarkLine = isDarkLine;
	this.doCorrectPosition = doCorrectPosition;
	this.doEstimateWidth = doEstimateWidth;
	this.doExtendLine = doExtendLine;
	junctions = new Junctions(ip.getSliceNumber());
	lines = get_lines(sigma, upperThresh, lowerThresh, minLength, maxLength, ip.getHeight(),
			ip.getWidth(), ip, junctions, overlapOption);
	return lines;
}
	
	
	private void assignLinesToJunctions(Lines lines, Junctions junctions){
		for (Junction j : junctions) {
			j.lineCont1 = lines.get(j.cont1);
			j.lineCont2 = lines.get(j.cont2);
		}
	}

	public Options getUsedParamters() {
		return opts;
	}

	public Junctions getJunctions() {
		return junctions;
	}
	
	private void addAdditionalJunctionPointsAndLines(Lines lines, Junctions junctions){
		
		for(int i = 0; i < junctions.size(); i++){
			Junction splitPoint = junctions.get(i); //Split point!
		
			log("Process Splitpoint " + splitPoint.getLine1().getID() +"-"+splitPoint.getLine2().getID() + " Pos: " + splitPoint.pos);
		//splitPoint.pos!=0&&splitPoint.pos!=(splitPoint.getLine1().num-1)
			if(!alreadyProcessedJunctionPoints.contains(i)){
			
				
				/*
				 * Find Junctions with the same position as the split point
				 */
				Junctions junctionsWithTheSamePosition = new Junctions(junctions.getFrame());
				alreadyProcessedJunctionPoints.add(i);
				junctionsWithTheSamePosition.add(splitPoint);
				for(int j = i+1; j < junctions.size(); j++){
					if(!alreadyProcessedJunctionPoints.contains(j)){
						Junction junc2 = junctions.get(j);
						if(Math.abs(junc2.x-splitPoint.x) < 0.01 && Math.abs(junc2.y-splitPoint.y) < 0.01){
							alreadyProcessedJunctionPoints.add(j);
							junctionsWithTheSamePosition.add(junc2);
						}
					}
				}
				
				/*
				 * Connect all lines which are connected with the processed line also with each other (new junctions point)
				 */
				
				ArrayList<Line> connectedWithProcessedLine = new ArrayList<Line>();
				ArrayList<Integer> connectedWithProcessedIndex = new ArrayList<Integer>();
				for (Junction junc : junctionsWithTheSamePosition) {
					connectedWithProcessedLine.add(junc.getLine2());
					connectedWithProcessedIndex.add(junc.cont2);
				}
				for(int j = 0; j < connectedWithProcessedLine.size(); j++){
					for(int k = j+1; k < connectedWithProcessedLine.size(); k++){
						Line l1 = connectedWithProcessedLine.get(j);
						Line l2 = connectedWithProcessedLine.get(k);
						Junction junc = new Junction();
						junc.lineCont1 = l1;
						junc.lineCont2 = l2;
						junc.x = splitPoint.x;
						junc.y = splitPoint.y;
						junc.cont1 = connectedWithProcessedIndex.get(j);
						junc.cont2 = connectedWithProcessedIndex.get(k);
						junc.pos = l1.getStartOrdEndPosition(junc.x, junc.y);
						junctions.add(junc); 
						log("Connect " + junc.getLine1().getID() +"-"+junc.getLine2().getID() + " Pos: " + junc.pos);
						//l1.setContourClass(reconstructContourClass(l1, l1.getStartOrdEndPosition(junc.x, junc.y)));
					//	l2.setContourClass(reconstructContourClass(l2, l2.getStartOrdEndPosition(junc.x, junc.y)));
						alreadyProcessedJunctionPoints.add(junctions.size()-1);
					}
				}
				
				/*
				 * Split the line in two line at the split point if it is not at the end or beginning of the line
				 */
				Line l1 = splitPoint.getLine1();
				int pos = splitPoint.pos;
				boolean isClosedContour = l1.col[0] == l1.col[l1.num-1] && l1.row[0] == l1.row[l1.num-1];
				
				if(isClosedContour){
					l1.setContourClass(LinesUtil.contour_class.cont_closed);
					l1.setContourClass(reconstructContourClass(l1, l1.getStartOrdEndPosition(splitPoint.x, splitPoint.y)));
				}
				log("Pos: " + pos + " num: " + l1.num);
				if(pos!=0 && pos != (l1.num-1) && !isClosedContour){
					//All data up to pos (included)
					int keepLength = pos+1;
					
				
					float[] keepAsymmetry = new float[keepLength];
					float[] keepIntensity = new float[keepLength];
						
					float[] keepAngle = new float[keepLength];
					float[] keepWidth_l = new float[keepLength];
					float[] keepWidth_r = new float[keepLength];
				
					
					float[] keepCol = new float[keepLength];
					float[] keepRow = new float[keepLength];
					float[] keepResponse = new float[keepLength];
					
					
					
					//All data from pos (included)
					int splitSize = (int) (l1.num-pos);
					
					float[] splitAsymmetry = new float[splitSize];
					float[] splitIntensity = new float[splitSize];
					
					float[] splitAngle = new float[splitSize];
					float[] splitWidth_l = new float[splitSize];
					float[] splitWidth_r = new float[splitSize];
					
					float[] splitCol = new float[splitSize];
					float[] splitRow = new float[splitSize];
					float[] splitResponse = new float[splitSize];
					
					
					//Copy data
					if(doEstimateWidth){
						if(doCorrectPosition){
							System.arraycopy(l1.asymmetry, 0, keepAsymmetry, 0, keepLength);
							System.arraycopy(l1.asymmetry, pos, splitAsymmetry, 0, splitSize);
						
							System.arraycopy(l1.intensity, 0, keepIntensity, 0, keepLength);
							System.arraycopy(l1.intensity, pos, splitIntensity, 0, splitSize);
						}
						
						System.arraycopy(l1.angle, 0, keepAngle, 0, keepLength);
						System.arraycopy(l1.angle, pos, splitAngle, 0, splitSize);
						
						System.arraycopy(l1.width_l, 0, keepWidth_l, 0, keepLength);
						System.arraycopy(l1.width_l, pos, splitWidth_l, 0, splitSize);
						
						System.arraycopy(l1.width_r, 0, keepWidth_r, 0, keepLength);
						System.arraycopy(l1.width_r, pos, splitWidth_r, 0, splitSize);
					}
					
					System.arraycopy(l1.col, 0, keepCol, 0, keepLength);
					System.arraycopy(l1.col, pos, splitCol, 0, splitSize);
					
					System.arraycopy(l1.row, 0, keepRow, 0, keepLength);
					System.arraycopy(l1.row, pos, splitRow, 0, splitSize);
					
					System.arraycopy(l1.response, 0, keepResponse, 0, keepLength);
					System.arraycopy(l1.response, pos, splitResponse, 0, splitSize);
					
					
					
					
					
					
					
					//Generate new line
					Line lNew = new Line();
					lNew.angle = splitAngle;
					lNew.asymmetry = splitAsymmetry;
					lNew.col = splitCol;
					lNew.row = splitRow;
					lNew.response = splitResponse;
					lNew.intensity = splitIntensity;
					lNew.width_l = splitWidth_l;
					lNew.width_r = splitWidth_r;
					lNew.num = splitSize;
					lNew.setContourClass(l1.getContourClass());
					lNew.setFrame(l1.getFrame());
					lines.add(lNew);
					int newID = lNew.getID();
					
					//Update junctions
					
					//Add additional junction points for the split point
					Set<Integer> lineIds = new HashSet<Integer>(); //All IDs which are connected at this junction point
					for (Junction junc : junctionsWithTheSamePosition) {
						lineIds.add(junc.getLine1().getID());
						lineIds.add(junc.getLine2().getID());
					}
					Iterator<Integer> idIt = lineIds.iterator();
					while (idIt.hasNext()) {
						int id = idIt.next();
						int connectWithLineID =  lines.getIndexByID(id);
						Line connectWith = lines.get(connectWithLineID);
						Junction j = new Junction();
						j.cont1 = lines.size()-1;
						j.cont2 = connectWithLineID;
						j.lineCont1 = lNew;
						j.lineCont2 = connectWith;
						j.x = splitPoint.x;
						j.y = splitPoint.y;
						j.pos = lNew.getStartOrdEndPosition(splitPoint.x, splitPoint.y);
						//lNew.setContourClass(reconstructContourClass(lNew, j.pos));
						//connectWith.setContourClass(reconstructContourClass(connectWith, connectWith.getStartOrdEndPosition(splitPoint.x, splitPoint.y)));
						junctions.add(j);
						log("Connect " + j.getLine1().getID() +"-"+j.getLine2().getID() + " Pos: " + j.pos);
						alreadyProcessedJunctionPoints.add(junctions.size()-1);
					}

					//Update following junctions point
					for(int j = 0; j < junctions.size(); j++){
						Junction junc2 = junctions.get(j);
						if(junc2.cont1 == splitPoint.cont1 && junc2.pos>splitPoint.pos){
							log("Update From " + junc2.getLine1().getID() +"-"+junc2.getLine2().getID() + " Pos: " + junc2.pos);
							junc2.cont1 = lines.getIndexByID(newID);
							junc2.lineCont1 = lNew;
							junc2.pos = junc2.pos-splitPoint.pos;
							log("Update To " + junc2.getLine1().getID() +"-"+junc2.getLine2().getID() + " Pos: " + junc2.pos);
						}
						
						double[] min = minDistance(junc2.getLine2(), junc2.x, junc2.y);
						
						if(junc2.cont2 == splitPoint.cont1 && ((int)min[1])>splitPoint.pos){
							junc2.cont2 = lines.getIndexByID(newID);
							junc2.lineCont2 = lNew;
						}
						
						
						
					}
					
					//Update Line 1
					//Overwrite line data
					l1.angle = keepAngle;
					l1.asymmetry = keepAsymmetry;
					l1.col = keepCol;
					l1.row = keepRow;
					l1.response = keepResponse;
					l1.intensity = keepIntensity;
					l1.width_l = keepWidth_l;
					l1.width_r = keepWidth_r;
					l1.num = keepLength;
					
					//Update position of splitpoint
					log("Set Splitpoint Position from " + splitPoint.pos);
					splitPoint.pos = l1.getStartOrdEndPosition(splitPoint.x, splitPoint.y);
					log("Set Splitpoint Position to " + splitPoint.pos);
					lines.set(splitPoint.cont1, l1);
					
				}
				
				
			}
		}
	}

	private Junctions fixJunctions(Lines lines, Junctions junctions) {
		/*
		 * For some reason, the x and y coordinates are permuted
		 */
		for (Junction junction : junctions) {
			float help = junction.x;
			junction.x = junction.y;
			junction.y = help;
		}

		Junctions newJunctions = new Junctions(junctions.getFrame());
		ArrayList<Point2D.Float> processedJunctions = new ArrayList<Point2D.Float>();
		for(int i = 0; i < junctions.size(); i++){
			Junction junc = junctions.get(i);
			Line mainLine = null;
			int mainLineIndex = -1;
			int mainLinePos = -1;
			ArrayList<Line> secondaryLines = new ArrayList<Line>();
			ArrayList<Integer> secondaryLineIndex = new ArrayList<Integer>();
			ArrayList<Integer> secondaryLinePos = new ArrayList<Integer>();
			
			//Verarbeite jede Junction-Position nur einmal.
			if(!processedJunctions.contains(new Point2D.Float(junc.x, junc.y))){ //processed[(int)junc.x][(int)junc.y]==0

				processedJunctions.add(new Point2D.Float(junc.x, junc.y));
				
				/*
				 * Finde die Sekundärlinien und Hauptlinien
				 */
				for(int j = 0; j < lines.size(); j++) {
					Line l = lines.get(j);
				
					double[] mindist = minDistance(l, junc.x, junc.y);
					if(mindist[0]<0.1){ //Wenn der Punkt auf der Linie liegt, analysiere genauer
						
						if(mindist[1]==0 || mindist[1]==(l.num-1)){ //Wenn der Junction-Point am Ende oder am Anfang liegt, ist es sekundäre Linie.
							secondaryLines.add(l);
							secondaryLineIndex.add(j);
							secondaryLinePos.add((int)mindist[1]);
						} else {			// Wenn er innerhalb der Linie liegt, ist dies die Hauptlinie.
			
							if(mainLine!=null){
								if(mainLine.getID()==l.getID()){
									continue;
								}
								log("Äh, zwei Hauptlininen geht nich..." + mainLine.getID() + " x " + junc.x + " y " + junc.y);
								log("Äh, zwei Hauptlininen geht nich..." + l.getID() + " x " + junc.x + " y " + junc.y);
							}
							mainLine = l;
							mainLineIndex = j;
							mainLinePos = (int) mindist[1];
							
							
						}
					}
				}
				if(mainLine!=null){
					for (int j = 0; j < secondaryLines.size(); j++) {
						Junction newJunc = new Junction();
						newJunc.cont1 = mainLineIndex;
						newJunc.cont2 = secondaryLineIndex.get(j);
						newJunc.x = junc.x;
						newJunc.y = junc.y;
						newJunc.pos = mainLinePos;
						//lines.get(newJunc.cont1).setContourClass(reconstructContourClass(lines.get(newJunc.cont1), mainLinePos));
						//lines.get(newJunc.cont2).setContourClass(reconstructContourClass(lines.get(newJunc.cont2), secondaryLinePos.get(j)));
						newJunctions.add(newJunc);
						log("NewJunc Mainline: " + lines.get(newJunc.cont1).getID() + "-" + lines.get(newJunc.cont2).getID() + " pos " + newJunc.pos + " num " + lines.get(newJunc.cont1).num);
						
					}
				}else{
					//In manchen Fällen gibt es keine Hauptlinie... (bug im Algorithmus, ich bin aber nicht fähig ihn zu finden=.
					HashSet<Integer> uniqueIDs = new HashSet<Integer>();
					ArrayList<Line> uniqueLines = new ArrayList<Line>();
					ArrayList<Integer> uniqueLineIndex = new ArrayList<Integer>();
					ArrayList<Integer> uniqueLinePos = new ArrayList<Integer>();
					for (int j = 0; j < secondaryLines.size(); j++) {
						if(!uniqueIDs.contains(secondaryLines.get(j).getID())){
							uniqueIDs.add(secondaryLines.get(j).getID());
							uniqueLines.add(secondaryLines.get(j));
							uniqueLineIndex.add(secondaryLineIndex.get(j));
							uniqueLinePos.add(secondaryLinePos.get(j));
							
						}
						
					}
					for(int j = 0; j < uniqueLines.size(); j++){
						for(int k = j+1; k < uniqueLines.size(); k++){
							Junction newJunc = new Junction();
							newJunc.cont1 = uniqueLineIndex.get(j);
							newJunc.cont2 = uniqueLineIndex.get(k);
							newJunc.x = junc.x;
							newJunc.y = junc.y;
							newJunc.pos = uniqueLinePos.get(j);
							newJunctions.add(newJunc);
							log("NewJunc Second: " + lines.get(newJunc.cont1).getID() + "-" + lines.get(newJunc.cont2).getID() + " pos " + newJunc.pos + " num " + lines.get(newJunc.cont1).num);

							//lines.get(newJunc.cont1).setContourClass(reconstructContourClass(lines.get(newJunc.cont1), uniqueLinePos.get(j)));
							//lines.get(newJunc.cont2).setContourClass(reconstructContourClass(lines.get(newJunc.cont2), uniqueLinePos.get(k)));
							alreadyProcessedJunctionPoints.add(newJunctions.size()-1);
					
							
						}
					}
				}
			
			}
		}
		return newJunctions;
		
	}
	
	private LinesUtil.contour_class reconstructContourClass(Line l, int pos){
		LinesUtil.contour_class currentClass = l.getLineClass();
	
		boolean hasJunctionAtStartpoint = pos==0?true:false;
		boolean hasJunctionAtEndpoint = pos==(l.num-1)?true:false;
		
		if(currentClass == LinesUtil.contour_class.cont_no_junc && hasJunctionAtStartpoint){
			return LinesUtil.contour_class.cont_start_junc;
		}
		
		if(currentClass == LinesUtil.contour_class.cont_no_junc && hasJunctionAtEndpoint){
			return LinesUtil.contour_class.cont_end_junc;
		}
		
		if(currentClass == LinesUtil.contour_class.cont_start_junc && hasJunctionAtEndpoint){
			return LinesUtil.contour_class.cont_both_junc;
		}
		
		if(currentClass == LinesUtil.contour_class.cont_end_junc && hasJunctionAtStartpoint){
			return LinesUtil.contour_class.cont_both_junc;
		}
		
		if(currentClass == LinesUtil.contour_class.cont_closed && (hasJunctionAtEndpoint || hasJunctionAtStartpoint)){
			return LinesUtil.contour_class.cont_both_junc;
		}
		
		return currentClass;
		
	}
	/**
	 * 
	 * @param l Line
	 * @param x x-Position
	 * @param y y-Position
	 * @return Double Array [0] minimal distance [1] position of minimal distance
	 */
	private double[] minDistance(Line l, float x, float y){
		double min = Double.MAX_VALUE;
		double index = -1;
		for(int i = 0; i < l.num; i++){
			double d = Math.sqrt(Math.pow(l.col[i]-x, 2)+Math.pow(l.row[i]-y, 2));
			if(d< min){
				min =d;
				index = i;
			}
		}
		double[] ret =  {min,index};
		return ret;
	}
	
/*	private void deleteContour(Lines contours, Junctions junctions, Line c) {

		ArrayList<Junction> remove = new ArrayList<Junction>();
		for (Junction junction : junctions) {

			if (contours.get((int) junction.cont1).getID() == c.getID()
					|| contours.get((int) junction.cont2).getID() == c.getID()) {
				remove.add(junction);
			}

		}
		for (Junction junction : remove) {
			junctions.remove(junction);
		}

		contours.remove(c);
	} */

	//To be removed once the problem with junctions.cont1 & .cont2 is solved.
	private void deleteJunctions(Lines contours, Junctions junctions, Line c) {
		deleteJunctions(contours, junctions, c, OverlapOption.NONE);
	}
	private void deleteJunctions(Lines contours, Junctions junctions, Line c, OverlapOption overlapOption) {	

		ArrayList<Junction> remove = new ArrayList<Junction>();
		for (Junction junction : junctions) {
			// This if() should be removed once cont1 and 2 contain the same info whatever OverlapOption
			if (overlapOption == OverlapOption.SLOPE) {
				if (junction.cont1 == c.getID() || junction.cont2 == c.getID())  {
					log("Removing junction between line IDs"+ junction.cont1+ " and "+junction.cont2);
					remove.add(junction);
				}
			} else {
				if (contours.get((int) junction.cont1).getID() == c.getID()
					|| contours.get((int) junction.cont2).getID() == c.getID())  {
					log("Removing junction between line idx "+ junction.cont1+ " and "+junction.cont2);
					remove.add(junction);
				}	
			}
		}
		for (Junction junction : remove) {
			junctions.remove(junction);
		}
	}

	private void fixContours(Lines contours, Junctions junctions) {

		ArrayList<Line> remove = new ArrayList<Line>();
		// Contours with only a single position cant be valid.
		for (Line contour : contours) {
			if (contour.num == 1) {
				deleteJunctions(contours,junctions,contour);
				remove.add(contour);
			}
			//If the results are corrupted, this informationen has to be reconstructed in fixJunctions
			contour.setContourClass(LinesUtil.contour_class.cont_no_junc);
		}

		for (Line c : remove) {
       		contours.remove(c);
       	}

		// For some reason the first and the last element are the same. Delete
		// it!
		
		if (contours.size() >= 2) {
			if(contours.get(0).getID() == contours.get(contours.size()-1).getID()){
				contours.remove(contours.size() - 1);	
			}
		}
	}
    
    private void pruneContours(Lines contours, Junctions junctions, double minLength, double maxLength, OverlapOption overlapOption) {
      ArrayList<Line> remove = new ArrayList<Line>();

      log("Pruning lines:");
      for (Line c : contours) {
      	if ((c.estimateLength() < minLength) || (maxLength > 0 && c.estimateLength() > maxLength)) {
        	log("Removing line "+c.getID()+ " of length "+c.estimateLength());
        	deleteJunctions(contours, junctions, c,overlapOption);
        	remove.add(c);
       	} else {
       		log("Keeping line "+c.getID()+ " of length "+c.estimateLength());
       	}
       }
       for (Line c : remove) {
       	
       	contours.remove(c);
       }
    }
    
    
    
    private Lines get_lines(double sigma, double high, double low, double minLength, double maxLength, int rows,
			int cols, ImageProcessor in_img, Junctions resultJunction, OverlapOption overlapOption) {
		FloatProcessor image;
		Lines contours = new Lines(in_img.getSliceNumber());
		int num_cont = 0;
		opts = new Options(-1.0, -1.0, -1.0, isDarkLine ? LinesUtil.MODE_DARK
				: LinesUtil.MODE_LIGHT, -1.0, -1.0, doCorrectPosition, doEstimateWidth,
				doExtendLine, false, false, false, overlapOption);

		opts.sigma = sigma;
		opts.high = high;
		opts.low = low;
		check_sigma(opts.sigma, cols, rows);

		opts.minLength = minLength; 
		opts.maxLength = maxLength; 

		OverlapResolver resolver = null;

		switch (overlapOption) {
			default:
			case NONE:
				break;
			case SLOPE: resolver = new SlopeOverlapResolver();
				break;
		}

		int i2, j2;
		// //(float *) malloc(rows*cols*sizeof(float));
		float[] imgpxls = new float[cols * rows];
		for (i2 = 0; i2 < rows; i2++)
			for (j2 = 0; j2 < cols; j2++)
				imgpxls[i2 * cols + j2] = in_img.getf(j2, i2);
		image = new FloatProcessor(cols, rows, imgpxls);
		MutableInt hnum_cont = new MutableInt(num_cont);
		float[] imgpxls2 = (float[]) image.getPixels();
		Position p = new Position();
		p.detect_lines(imgpxls2, cols, rows, contours, hnum_cont, opts.sigma,
				opts.low, opts.high, opts.mode, opts.width, opts.correct,
				opts.extend, resultJunction);
		num_cont = hnum_cont.getValue();

	//	lines = contours;
		fixContours(contours,resultJunction);
		alreadyProcessedJunctionPoints = new HashSet<Integer>();
		//Reconstruct solution from junction points. This have to be done, because in raw cases
		//the algorithm corrupts the results. However, I was not able to find that bug so I decided
		//to reconstruct the solution from the information which were not be corrupted.

		resultJunction = fixJunctions(contours,resultJunction);

		assignLinesToJunctions(contours,resultJunction);
		
		addAdditionalJunctionPointsAndLines(contours,resultJunction);
		Collections.sort(resultJunction);
		junctions = resultJunction;

		/*
		 * RECONSTRUCTION OF CONTOUR CLASS
		 */
		//Reset contour class
		for(int i = 0; i < contours.size(); i++){
			contours.get(i).setContourClass(LinesUtil.contour_class.cont_no_junc);
		}
		
		//Find closed lines
		for(int i = 0; i < contours.size(); i++){
			boolean isClosedContour = contours.get(i).col[0] == contours.get(i).col[contours.get(i).num-1] && contours.get(i).row[0] == contours.get(i).row[contours.get(i).num-1];
			if(isClosedContour){
				contours.get(i).setContourClass(LinesUtil.contour_class.cont_closed);	
			}
		}

	    //Reconstruction contour class
		for(int i = 0; i < junctions.size(); i++){
			Junction j = junctions.get(i);
			j.getLine1().setContourClass(reconstructContourClass(j.getLine1(), j.pos));
			float x = j.getLine1().getXCoordinates()[j.pos];
			float y = j.getLine1().getYCoordinates()[j.pos];
			j.getLine2().setContourClass(reconstructContourClass(j.getLine2(),j.getLine2().getStartOrdEndPosition(x, y)));
		}
		
		if (resolver != null) contours = resolver.resolve(contours, junctions, bechatty);

		if (minLength != 0 || maxLength != 0) {
			pruneContours(contours, junctions, minLength, maxLength,overlapOption);
		}

		
		return contours;

	}
	private void log(String s){
		if(bechatty){
			IJ.log(s);
		}
	}
	private void check_sigma(double sigma, int width, int height) {
		int min_dim;
		min_dim = width < height ? width : height;
		if (sigma < 0.4)
			IJ.error(LinesUtil.ERR_SOR, "< 0.4");
		if (LinesUtil.MASK_SIZE(LinesUtil.MAX_SIZE_MASK_2, sigma) >= min_dim)
			IJ.error(LinesUtil.ERR_SOR, "too large for image size");
	}

}
