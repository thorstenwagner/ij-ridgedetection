package de.biomedical_imaging.ij.steger;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.commons.lang3.mutable.MutableLong;

import com.sun.tools.internal.xjc.reader.gbind.ConnectedComponent;

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
	private ImageProcessor ip;
	Set<Integer> alreadyProcessedJunctionPoints;

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
			double upperThresh, double lowerThresh, boolean isDarkLine,
			boolean doCorrectPosition, boolean doEstimateWidth,
			boolean doExtendLine) {
		this.ip = ip;
		this.isDarkLine = isDarkLine;
		this.doCorrectPosition = doCorrectPosition;
		this.doEstimateWidth = doEstimateWidth;
		this.doExtendLine = doExtendLine;
		junctions = new Junctions(ip.getSliceNumber());
		lines = get_lines(sigma, upperThresh, lowerThresh, ip.getHeight(),
				ip.getWidth(), ip, junctions);
		fixContours();
		alreadyProcessedJunctionPoints = new HashSet<Integer>();
		
		//Reconstruct solution from junction points. This have to be done, because in raw cases
		//the algorithm corrups the results. However, I was not able to find that bug so I decided
		//to reconstruct the solution from the information which were not corrupted.
		fixJunctions();
		assignLinesToJunctions();
		addAdditionalJunctionPointsAndLines();
		Collections.sort(junctions);
		return lines;
	}
	
	private void assignLinesToJunctions(){
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
	
	private void addAdditionalJunctionPointsAndLines(){
		
		for(int i = 0; i < junctions.size(); i++){
			Junction splitPoint = junctions.get(i); //Split point!
			if(!alreadyProcessedJunctionPoints.contains(i)){
				if(splitPoint.pos>=lines.get(splitPoint.cont1).num){
					alreadyProcessedJunctionPoints.add(i);
					continue;
				}
			
				
				/*
				 * Find Junctions with the same position
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
						junc.pos = l1.getIndexOfPosition(junc.x, junc.y);
						junctions.add(junc); 
						alreadyProcessedJunctionPoints.add(junctions.size()-1);
					}
				}
				
				/*
				 * Split the line in two line at the split point 
				 */
				Line l1 = splitPoint.getLine1();
				int pos = splitPoint.pos;
				
				if(pos!=0 && pos != (l1.num-1)){
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
					lNew.cont_class = l1.cont_class;
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
						j.pos = lNew.getIndexOfPosition(splitPoint.x, splitPoint.y);
						junctions.add(j);
						alreadyProcessedJunctionPoints.add(junctions.size()-1);
					}

					//Update following junctions point
					for(int j = 0; j < junctions.size(); j++){
						Junction junc2 = junctions.get(j);
						if(junc2.cont1 == splitPoint.cont1 && junc2.pos>splitPoint.pos){
							junc2.cont1 = lines.getIndexByID(newID);
							junc2.lineCont1 = lNew;
							junc2.pos = junc2.pos-splitPoint.pos;
						}
					}
					
					//Update position of splitpoint
					splitPoint.pos = l1.getIndexOfPosition(splitPoint.x, splitPoint.y);
					lines.set(splitPoint.cont1, l1);
					
					
				}
				
				
			}
		}
	}

	private void fixJunctions() {
		/*
		 * For some reason, the x and y coordinates are permuted
		 */
		for (Junction junction : junctions) {
			float help = junction.x;
			junction.x = junction.y;
			junction.y = help;
		}
		
		Junctions newJunctions = new Junctions(junctions.getFrame());
		int[][] processed = new int[ip.getWidth()][ip.getHeight()];
		for(int i = 0; i < junctions.size(); i++){
			Junction junc = junctions.get(i);
			Line mainLine = null;
			int mainLineIndex = -1;
			int mainLinePos = -1;
			ArrayList<Line> secondaryLines = new ArrayList<Line>();
			ArrayList<Integer> secondaryLineIndex = new ArrayList<Integer>();
			ArrayList<Integer> secondaryLinePos = new ArrayList<Integer>();
			if(processed[(int)junc.x][(int)junc.y]==0){
				processed[(int)junc.x][(int)junc.y]=1;
				for(int j = 0; j < lines.size(); j++) {
					Line l = lines.get(j);
				
					double[] mindist = minDistance(l, junc.x, junc.y);
					if(mindist[0]==0){
						if(mindist[1]==0 || mindist[1]==(l.num-1)){
							secondaryLines.add(l);
							secondaryLineIndex.add(j);
							secondaryLinePos.add((int)mindist[1]);
						}else {
			
							if(mainLine!=null){
								if(mainLine.getID()==l.getID()){
									continue;
								}
								IJ.error("Äh, zwei Hauptlininen geht nich..." + mainLine.getID() + " x " + junc.x + " y " + junc.y);
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
						newJunctions.add(newJunc);
						
					}
				}else{
					//In manchen Fällen gibt es keine Hauptlinie... (bug im Algorithmus, ich bin aber nicht fähig ihn zu finden.
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
							newJunc.cont1 = secondaryLineIndex.get(j);
							newJunc.cont2 = secondaryLineIndex.get(k);
							newJunc.x = junc.x;
							newJunc.y = junc.y;
							newJunc.pos = uniqueLinePos.get(j);
							newJunctions.add(newJunc);
							alreadyProcessedJunctionPoints.add(newJunctions.size()-1);
						}
					}
				}
			
			}
		}
		junctions = newJunctions;
		
	}

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
	private void deleteContour(Line c) {

		ArrayList<Junction> remove = new ArrayList<Junction>();
		for (Junction junction : junctions) {

			if (lines.get((int) junction.cont1).getID() == c.getID()
					|| lines.get((int) junction.cont2).getID() == c.getID()) {
				remove.add(junction);
			}

		}
		for (Junction junction : remove) {
			junctions.remove(junction);
		}

		lines.remove(c);
	}

	private void fixContours() {

		// Contours with only a single position cant be valid.
		for (Line contour : lines) {
			if (contour.num == 1) {
				deleteContour(contour);

			}
		}

		// For some reason the first and the last element are the same. Delete
		// it!
		if (lines.size() >= 2) {
			if(lines.get(0).getID() == lines.get(lines.size()-1).getID()){
				lines.remove(lines.size() - 1);	
			}
		}
		
		

	}

	private Lines get_lines(double sigma, double high, double low, int rows,
			int cols, ImageProcessor in_img, Junctions resultJunction) {
		FloatProcessor image;
		Lines contours = new Lines(in_img.getSliceNumber());
		int num_cont = 0;
		opts = new Options(-1.0, -1.0, -1.0, isDarkLine ? LinesUtil.MODE_DARK
				: LinesUtil.MODE_LIGHT, doCorrectPosition, doEstimateWidth,
				doExtendLine, false, false, false);

		opts.sigma = sigma;
		opts.high = high;
		opts.low = low;
		check_sigma(opts.sigma, cols, rows);

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

		return contours;

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
