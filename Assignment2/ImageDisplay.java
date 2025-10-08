
import java.awt.*;
import java.awt.image.*;
import java.io.*;
import javax.swing.*;


public class ImageDisplay {

	JFrame frame;
	JLabel lbIm1, lbIm2;
	BufferedImage imgOne, imgTwo;

	// Modify the height and width values here to read and display an image with
  	// different dimensions. 
	int width = 352;
	int height = 288;
	
	// DCT and quantization related variables
	double[][][] dctCoefficients;
	int quantizationLevel;
	int deliveryMode;
	int latency;
	
	// Progress tracking variables
	int totalBlocks;
	int totalCoefficients;
	int totalBitLevels;
	int currentProgress;

	/** Read Image RGB
	 *  Reads the image of given width and height at the given imgPath into the provided BufferedImage.
	 */
	private void readImageRGB(int width, int height, String imgPath, BufferedImage img)
	{
		try
		{
			int frameLength = width*height*3;

			File file = new File(imgPath);

            RandomAccessFile raf = new RandomAccessFile(file, "r");

			raf.seek(0);

			long len = frameLength;
			byte[] bytes = new byte[(int) len];

			raf.read(bytes);

			int ind = 0;
			for(int y = 0; y < height; y++)
			{
				for(int x = 0; x < width; x++)
				{
					byte a = 0;
					byte r = bytes[ind];
					byte g = bytes[ind+height*width];
					byte b = bytes[ind+height*width*2]; 

					int pix = 0xff000000 | ((r & 0xff) << 16) | ((g & 0xff) << 8) | (b & 0xff);
					//int pix = ((a << 24) + (r << 16) + (g << 8) + b);
					img.setRGB(x,y,pix);
					ind++;
				}
			}
		}
		catch (FileNotFoundException e) 
		{
			e.printStackTrace();
		} 
		catch (IOException e) 
		{
			e.printStackTrace();
		}
	}

	/**
	 * Apply DCT to an 8x8 block
	 */
	private double[][] applyDCT(double[][] block) {
		double[][] dct = new double[8][8];
		double alpha, beta;
		
		for (int u = 0; u < 8; u++) {
			for (int v = 0; v < 8; v++) {
				alpha = (u == 0) ? 1.0 / Math.sqrt(2) : 1.0;
				beta = (v == 0) ? 1.0 / Math.sqrt(2) : 1.0;
				
				double sum = 0.0;
				for (int x = 0; x < 8; x++) {
					for (int y = 0; y < 8; y++) {
						sum += block[x][y] * 
							   Math.cos(((2 * x + 1) * u * Math.PI) / 16) *
							   Math.cos(((2 * y + 1) * v * Math.PI) / 16);
					}
				}
				dct[u][v] = 0.25 * alpha * beta * sum;
			}
		}
		return dct;
	}

	/**
	 * Apply Inverse DCT to an 8x8 block
	 */
	private double[][] applyIDCT(double[][] dct) {
		double[][] block = new double[8][8];
		double alpha, beta;
		
		for (int x = 0; x < 8; x++) {
			for (int y = 0; y < 8; y++) {
				double sum = 0.0;
				for (int u = 0; u < 8; u++) {
					for (int v = 0; v < 8; v++) {
						alpha = (u == 0) ? 1.0 / Math.sqrt(2) : 1.0;
						beta = (v == 0) ? 1.0 / Math.sqrt(2) : 1.0;
						
						sum += alpha * beta * dct[u][v] *
							   Math.cos(((2 * x + 1) * u * Math.PI) / 16) *
							   Math.cos(((2 * y + 1) * v * Math.PI) / 16);
					}
				}
				block[x][y] = 0.25 * sum;
			}
		}
		return block;
	}

	/**
	 * Quantize DCT coefficients
	 */
	private double[][] quantize(double[][] dct) {
		double[][] quantized = new double[8][8];
		double quantStep = Math.pow(2, quantizationLevel);
		
		for (int u = 0; u < 8; u++) {
			for (int v = 0; v < 8; v++) {
				quantized[u][v] = Math.round(dct[u][v] / quantStep);
			}
		}
		return quantized;
	}

	/**
	 * Dequantize DCT coefficients
	 */
	private double[][] dequantize(double[][] quantized) {
		double[][] dct = new double[8][8];
		double quantStep = Math.pow(2, quantizationLevel);
		
		for (int u = 0; u < 8; u++) {
			for (int v = 0; v < 8; v++) {
				dct[u][v] = quantized[u][v] * quantStep;
			}
		}
		return dct;
	}

	/**
	 * Extract 8x8 block from image channel
	 */
	private double[][] extractBlock(int[][][] image, int channel, int blockX, int blockY) {
		double[][] block = new double[8][8];
		for (int x = 0; x < 8; x++) {
			for (int y = 0; y < 8; y++) {
				int pixelX = blockX * 8 + x;
				int pixelY = blockY * 8 + y;
				if (pixelX < width && pixelY < height) {
					block[x][y] = image[pixelY][pixelX][channel];
				}
			}
		}
		return block;
	}

	/**
	 * Place 8x8 block back into image channel
	 */
	private void placeBlock(double[][] block, int[][][] image, int channel, int blockX, int blockY) {
		for (int x = 0; x < 8; x++) {
			for (int y = 0; y < 8; y++) {
				int pixelX = blockX * 8 + x;
				int pixelY = blockY * 8 + y;
				if (pixelX < width && pixelY < height) {
					image[pixelY][pixelX][channel] = (int) Math.max(0, Math.min(255, block[x][y]));
				}
			}
		}
	}

	/**
	 * Convert BufferedImage to 3D array [y][x][channel]
	 */
	private int[][][] imageToArray(BufferedImage img) {
		int[][][] image = new int[height][width][3];
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				int rgb = img.getRGB(x, y);
				image[y][x][0] = (rgb >> 16) & 0xFF; // R
				image[y][x][1] = (rgb >> 8) & 0xFF;  // G
				image[y][x][2] = rgb & 0xFF;         // B
			}
		}
		return image;
	}

	/**
	 * Convert 3D array to BufferedImage
	 */
	private BufferedImage arrayToImage(int[][][] image) {
		BufferedImage img = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				int r = Math.max(0, Math.min(255, image[y][x][0]));
				int g = Math.max(0, Math.min(255, image[y][x][1]));
				int b = Math.max(0, Math.min(255, image[y][x][2]));
				int rgb = (r << 16) | (g << 8) | b;
				img.setRGB(x, y, rgb);
			}
		}
		return img;
	}

	/**
	 * Encode image using DCT and quantization
	 */
	private void encodeImage(BufferedImage original) {
		int[][][] image = imageToArray(original);
		int blocksX = (width + 7) / 8;
		int blocksY = (height + 7) / 8;
		
		// Initialize DCT coefficients storage
		dctCoefficients = new double[3][blocksY][blocksX * 64]; // 3 channels, each block has 64 coefficients
		
		// Process each channel
		for (int channel = 0; channel < 3; channel++) {
			for (int blockY = 0; blockY < blocksY; blockY++) {
				for (int blockX = 0; blockX < blocksX; blockX++) {
					// Extract 8x8 block
					double[][] block = extractBlock(image, channel, blockX, blockY);
					
					// Apply DCT
					double[][] dct = applyDCT(block);
					
					// Quantize
					double[][] quantized = quantize(dct);
					
					// Store coefficients in zigzag order
					int coeffIndex = blockX * 64;
					for (int u = 0; u < 8; u++) {
						for (int v = 0; v < 8; v++) {
							dctCoefficients[channel][blockY][coeffIndex + u * 8 + v] = quantized[u][v];
						}
					}
				}
			}
		}
	}

	/**
	 * Decode image based on delivery mode
	 */
	private void decodeImage() {
		int blocksX = (width + 7) / 8;
		int blocksY = (height + 7) / 8;
		int[][][] decodedImage = new int[height][width][3];
		
		// Initialize progress tracking
		totalBlocks = blocksX * blocksY;
		totalCoefficients = 64; // 8x8 DCT coefficients
		totalBitLevels = 8; // 1 to 8 bits
		currentProgress = 0;
		
		switch (deliveryMode) {
			case 1: // Baseline mode
				decodeBaseline(decodedImage, blocksX, blocksY);
				break;
			case 2: // Progressive spectral selection
				decodeProgressiveSpectral(decodedImage, blocksX, blocksY);
				break;
			case 3: // Progressive successive bit approximation
				decodeProgressiveBitApproximation(decodedImage, blocksX, blocksY);
				break;
		}
		
		imgTwo = arrayToImage(decodedImage);
		updateDisplay();
	}

	/**
	 * Baseline decoding - sequential block processing
	 */
	private void decodeBaseline(int[][][] decodedImage, int blocksX, int blocksY) {
		currentProgress = 0;
		
		for (int blockY = 0; blockY < blocksY; blockY++) {
			for (int blockX = 0; blockX < blocksX; blockX++) {
				// Process each channel
				for (int channel = 0; channel < 3; channel++) {
					// Get quantized coefficients
					double[][] quantized = new double[8][8];
					int coeffIndex = blockX * 64;
					for (int u = 0; u < 8; u++) {
						for (int v = 0; v < 8; v++) {
							quantized[u][v] = dctCoefficients[channel][blockY][coeffIndex + u * 8 + v];
						}
					}
					
					// Dequantize
					double[][] dct = dequantize(quantized);
					
					// Apply IDCT
					double[][] block = applyIDCT(dct);
					
					// Place block back in image
					placeBlock(block, decodedImage, channel, blockX, blockY);
				}
				
				// Update progress
				currentProgress++;
				
				// Update display after each block to show progressive reconstruction
				imgTwo = arrayToImage(decodedImage);
				updateDisplay();
				
				// Add latency to make the reconstruction process visible
				if (latency > 0) {
					try {
						Thread.sleep(latency);
					} catch (InterruptedException e) {
						e.printStackTrace();
					}
				}
			}
		}
	}

	/**
	 * Progressive spectral selection decoding
	 */
	private void decodeProgressiveSpectral(int[][][] decodedImage, int blocksX, int blocksY) {
		currentProgress = 0;
		
		// Start with DC coefficients only
		for (int coeff = 0; coeff < 64; coeff++) {
			for (int blockY = 0; blockY < blocksY; blockY++) {
				for (int blockX = 0; blockX < blocksX; blockX++) {
					// Process each channel
					for (int channel = 0; channel < 3; channel++) {
						// Get quantized coefficients up to current coefficient
						double[][] quantized = new double[8][8];
						int coeffIndex = blockX * 64;
						
						// Set coefficients based on current progressive level
						for (int u = 0; u < 8; u++) {
							for (int v = 0; v < 8; v++) {
								int currentCoeff = u * 8 + v;
								if (currentCoeff <= coeff) {
									quantized[u][v] = dctCoefficients[channel][blockY][coeffIndex + currentCoeff];
								} else {
									quantized[u][v] = 0; // Set unused coefficients to 0
								}
							}
						}
						
						// Dequantize
						double[][] dct = dequantize(quantized);
						
						// Apply IDCT
						double[][] block = applyIDCT(dct);
						
						// Place block back in image
						placeBlock(block, decodedImage, channel, blockX, blockY);
					}
				}
			}
			
			// Update progress
			currentProgress = coeff + 1;
			
			// Update display after each coefficient level to show progressive reconstruction
			imgTwo = arrayToImage(decodedImage);
			updateDisplay();
			
			// Add latency to make the reconstruction process visible
			if (latency > 0) {
				try {
					Thread.sleep(latency);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
		}
	}

	/**
	 * Progressive successive bit approximation decoding
	 */
	private void decodeProgressiveBitApproximation(int[][][] decodedImage, int blocksX, int blocksY) {
		currentProgress = 0;
		
		// Start with 1 bit, then 2 bits, etc.
		for (int bits = 1; bits <= 8; bits++) {
			for (int blockY = 0; blockY < blocksY; blockY++) {
				for (int blockX = 0; blockX < blocksX; blockX++) {
					// Process each channel
					for (int channel = 0; channel < 3; channel++) {
						// Get quantized coefficients
						double[][] quantized = new double[8][8];
						int coeffIndex = blockX * 64;
						for (int u = 0; u < 8; u++) {
							for (int v = 0; v < 8; v++) {
								double coeff = dctCoefficients[channel][blockY][coeffIndex + u * 8 + v];
								// Apply bit approximation
								quantized[u][v] = applyBitApproximation(coeff, bits);
							}
						}
						
						// Dequantize
						double[][] dct = dequantize(quantized);
						
						// Apply IDCT
						double[][] block = applyIDCT(dct);
						
						// Place block back in image
						placeBlock(block, decodedImage, channel, blockX, blockY);
					}
				}
			}
			
			// Update progress
			currentProgress = bits;
			
			// Update display after each bit level to show progressive reconstruction
			imgTwo = arrayToImage(decodedImage);
			updateDisplay();
			
			// Add latency to make the reconstruction process visible
			if (latency > 0) {
				try {
					Thread.sleep(latency);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
		}
	}

	/**
	 * Apply bit approximation to coefficient
	 */
	private double applyBitApproximation(double coeff, int bits) {
		if (coeff == 0) return 0;
		
		// Get sign
		double sign = coeff >= 0 ? 1 : -1;
		coeff = Math.abs(coeff);
		
		// Find most significant bit
		int msb = 0;
		double temp = coeff;
		while (temp >= 1) {
			temp /= 2;
			msb++;
		}
		
		// Apply bit approximation
		double result = 0;
		for (int i = 0; i < bits && (msb - i) >= 0; i++) {
			if (coeff >= Math.pow(2, msb - i)) {
				result += Math.pow(2, msb - i);
				coeff -= Math.pow(2, msb - i);
			}
		}
		
		return sign * result;
	}

	/**
	 * Update display with current decoded image and progress
	 */
	private void updateDisplay() {
		if (imgTwo != null) {
			lbIm2.setIcon(new ImageIcon(imgTwo));
			
			// Update progress text based on delivery mode
			String progressText = getProgressText();
			lbIm2.setText("Decoded Image - " + progressText);
			lbIm2.setHorizontalTextPosition(JLabel.CENTER);
			lbIm2.setVerticalTextPosition(JLabel.BOTTOM);
			
			frame.repaint();
		}
	}
	
	/**
	 * Get progress text based on current delivery mode and progress
	 */
	private String getProgressText() {
		double percentage = 0.0;
		String modeText = "";
		
		switch (deliveryMode) {
			case 1: // Baseline mode
				percentage = (double) currentProgress / totalBlocks * 100.0;
				modeText = String.format("Blocks: %d/%d (%.1f%%)", currentProgress, totalBlocks, percentage);
				break;
			case 2: // Progressive spectral selection
				percentage = (double) currentProgress / totalCoefficients * 100.0;
				modeText = String.format("Coefficients: %d/%d (%.1f%%)", currentProgress, totalCoefficients, percentage);
				break;
			case 3: // Progressive bit approximation
				percentage = (double) currentProgress / totalBitLevels * 100.0;
				modeText = String.format("Bit Level: %d/%d (%.1f%%)", currentProgress, totalBitLevels, percentage);
				break;
		}
		
		return modeText;
	}

	public void showIms(String[] args){

		// Parse command line arguments
		if (args.length != 4) {
			System.out.println("Usage: java ImageDisplay <inputImage> <quantizationLevel> <deliveryMode> <latency>");
			System.out.println("  quantizationLevel: 0-7 (compression level)");
			System.out.println("  deliveryMode: 1=baseline, 2=progressive spectral, 3=progressive bit approximation");
			System.out.println("  latency: milliseconds delay to make decoding process visible (0=instant)");
			System.exit(1);
		}

		String inputImage = args[0];
		quantizationLevel = Integer.parseInt(args[1]);
		deliveryMode = Integer.parseInt(args[2]);
		latency = Integer.parseInt(args[3]);

		// Validate parameters
		if (quantizationLevel < 0 || quantizationLevel > 7) {
			System.err.println("Quantization level must be between 0 and 7");
			System.exit(1);
		}
		if (deliveryMode < 1 || deliveryMode > 3) {
			System.err.println("Delivery mode must be 1, 2, or 3");
			System.exit(1);
		}
		if (latency < 0) {
			System.err.println("Latency must be non-negative");
			System.exit(1);
		}

		// Read in the specified image
		imgOne = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		readImageRGB(width, height, inputImage, imgOne);

		// Create frame with dual display
		frame = new JFrame("DCT Coder-Decoder - Original (Left) | Decoded (Right)");
		GridBagLayout gLayout = new GridBagLayout();
		frame.getContentPane().setLayout(gLayout);

		// Original image label
		lbIm1 = new JLabel(new ImageIcon(imgOne));
		lbIm1.setText("Original Image");
		lbIm1.setHorizontalTextPosition(JLabel.CENTER);
		lbIm1.setVerticalTextPosition(JLabel.BOTTOM);

		// Decoded image label (initially empty)
		lbIm2 = new JLabel("Decoded Image");
		lbIm2.setHorizontalTextPosition(JLabel.CENTER);
		lbIm2.setVerticalTextPosition(JLabel.BOTTOM);

		// Add labels to frame
		GridBagConstraints c = new GridBagConstraints();
		c.fill = GridBagConstraints.HORIZONTAL;
		c.anchor = GridBagConstraints.CENTER;
		c.weightx = 0.5;
		c.gridx = 0;
		c.gridy = 0;
		frame.getContentPane().add(lbIm1, c);

		c.gridx = 1;
		c.gridy = 0;
		frame.getContentPane().add(lbIm2, c);

		frame.pack();
		frame.setVisible(true);

		// Start encoding and decoding process
		System.out.println("Starting encoding process...");
		encodeImage(imgOne);
		System.out.println("Encoding complete. Starting decoding process...");
		
		// Start decoding in a separate thread to avoid blocking the UI
		Thread decodeThread = new Thread(() -> {
			decodeImage();
			System.out.println("Decoding complete.");
		});
		decodeThread.start();
	}

	public static void main(String[] args) {
		ImageDisplay ren = new ImageDisplay();
		ren.showIms(args);
	}

}
