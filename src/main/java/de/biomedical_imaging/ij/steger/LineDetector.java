package de.biomedical_imaging.ij.steger;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

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
		this.isDarkLine = isDarkLine;
		this.doCorrectPosition = doCorrectPosition;
		this.doEstimateWidth = doEstimateWidth;
		this.doExtendLine = doExtendLine;
		junctions = new Junctions(ip.getSliceNumber());
		lines = get_lines(sigma, upperThresh, lowerThresh, ip.getHeight(),
				ip.getWidth(), ip, junctions);
		fixContours();
		fixJunctions();
		assignLinesToJunctions();
		addAdditionalJunctionPointsAndLines();
		Collections.sort(junctions);
		return lines;
	}
	
	private void assignLinesToJunctions(){
		for (Junction j : junctions) {
			j.lineCont1 = lines.get((int) j.cont1);
			j.lineCont2 = lines.get((int) j.cont2);
		}
	}

	public Options getUsedParamters() {
		return opts;
	}

	public Junctions getJunctions() {
		return junctions;
	}
	
	private void addAdditionalJunctionPointsAndLines(){
		Set<Integer> alreadyProcessed = new HashSet<Integer>();
		for(int i = 0; i < junctions.size(); i++){
			Junction splitPoint = junctions.get(i); //Split point!
			if(!alreadyProcessed.contains(i)){
				
				/*
				 * Find Junctions with the same position
				 */
				Junctions junctionsWithTheSamePosition = new Junctions(junctions.getFrame());
				alreadyProcessed.add(i);
				junctionsWithTheSamePosition.add(splitPoint);
				for(int j = i+1; j < junctions.size(); j++){
					if(!alreadyProcessed.contains(j)){
						Junction junc2 = junctions.get(j);
						if(Math.abs(junc2.x-splitPoint.x) < 0.01 && Math.abs(junc2.y-splitPoint.y) < 0.01){
							alreadyProcessed.add(j);
							junctionsWithTheSamePosition.add(junc2);
						}
					}
				}
				
				/*
				 * Connect all lines which are connected with the processed line also with each other (new junctions point)
				 */
				ArrayList<Line> connectedWithProcessedLine = new ArrayList<Line>();
				ArrayList<Long> connectedWithProcessedIndex = new ArrayList<Long>();
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
						junctions.add(junc); 
						alreadyProcessed.add(junctions.size()-1);
					}
				}
				
				/*
				 * Go along the processed line and cut the line in two parts at every junction point on that line.
				 */
				Line l1 = splitPoint.getLine1();
				int pos = (int) splitPoint.pos;
				
				if(pos!=0 && pos != (l1.getNumber()-1)){
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
						j.cont1 = connectWithLineID;
						j.cont2 = lines.size()-1;
						j.lineCont1 = connectWith;
						j.lineCont2 = lNew;
						j.x = splitPoint.x;
						j.y = splitPoint.y;
						junctions.add(j);
						alreadyProcessed.add(junctions.size()-1);
					}

					//Update following junctions point
					for(int j = 0; j < junctions.size(); j++){
						Junction junc2 = junctions.get(j);
						if(junc2.cont1 == splitPoint.cont1 && junc2.pos>splitPoint.pos){
							//Junctionspoint which appear later on that line
							if(junc2.cont1 == splitPoint.cont1){
								junc2.cont1 = lines.getIndexByID(newID);
								junc2.lineCont1 = lNew;
								junc2.pos = junc2.pos-splitPoint.pos;
							}
						}
					}
					
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
			lines.remove(lines.size() - 1);
		}

	}

	private Lines get_lines(double sigma, double high, double low, int rows,
			int cols, ImageProcessor in_img, Junctions resultJunction) {
		FloatProcessor image;
		Lines contours = new Lines(in_img.getSliceNumber());
		long num_cont = 0;
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
		MutableLong hnum_cont = new MutableLong(num_cont);
		float[] imgpxls2 = (float[]) image.getPixels();
		Position p = new Position();
		p.detect_lines(imgpxls2, cols, rows, contours, hnum_cont, opts.sigma,
				opts.low, opts.high, opts.mode, opts.width, opts.correct,
				opts.extend, resultJunction);
		num_cont = hnum_cont.getValue();

		return contours;

	}

	private void check_sigma(double sigma, long width, long height) {
		long min_dim;
		min_dim = width < height ? width : height;
		if (sigma < 0.4)
			IJ.error(LinesUtil.ERR_SOR, "< 0.4");
		if (LinesUtil.MASK_SIZE(LinesUtil.MAX_SIZE_MASK_2, sigma) >= min_dim)
			IJ.error(LinesUtil.ERR_SOR, "too large for image size");
	}

}
