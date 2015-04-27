package de.biomedical_imaging.ij.steger;

import java.util.ArrayList;

import org.apache.commons.lang3.mutable.MutableLong;

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
		return lines;
	}

	public Options getUsedParamters() {
		return opts;
	}

	public Junctions getJunctions() {
		return junctions;
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
