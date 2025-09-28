#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <wx/dcbuffer.h>
#include <wx/wx.h>
#define DEBUG 0
using namespace std;
namespace fs = std::filesystem;

/**
 * Display an image using WxWidgets.
 * https://www.wxwidgets.org/
 */

/** Declarations*/

/**
 * Class that implements wxApp
 */
class MyApp : public wxApp {
public:
  bool OnInit() override;
};

/**
 * Class that implements wxFrame.
 * This frame serves as the top level window for the program
 */
class MyFrame : public wxFrame {
public:
  MyFrame(const wxString &title, string imagePath, int C, int M, int Q1, int Q2,
          int Q3);

private:
  void OnPaint(wxPaintEvent &event);
  wxImage inImage;
  wxImage outImage;
  wxScrolledWindow *scrolledWindow;
  int width;
  int height;
  int C, M, Q1, Q2, Q3;
};

/** Utility function to read image data */
unsigned char *readImageData(string imagePath, int width, int height);
double *from255to1(unsigned char *imageData, int width, int height);
unsigned char *from1to255(double *imageData, int width, int height);
unsigned char *quantize(unsigned char *imageData, int width, int height, int C,
                        int M, int Q1, int Q2, int Q3);
double *RGB2YUV(unsigned char *imageData, int width, int height);
unsigned char *YUV2RGB(double *imageData, int width, int height);
void saveToPngFile(string filePath, unsigned char *data, int width, int height);
void saveToRgbFile(string filePath, unsigned char *data, int width, int height);
double *nonUniformQuantizeYUV(const double *imageData, int width, int height,
                              int Q1, int Q2, int Q3, double Ymin = 0.0,
                              double Ymax = 1.0, double Umin = -0.436,
                              double Umax = 0.436, double Vmin = -0.615,
                              double Vmax = 0.615, int H = 256);
void testCases();
/** Definitions */

/**
 * Init method for the app.
 * Here we process the command line arguments and
 * instantiate the frame.
 */
bool MyApp::OnInit() {
  wxInitAllImageHandlers();

  // deal with command line arguments here
  cout << "Number of command line arguments: " << wxApp::argc << endl;
  if (wxApp::argc != 7) {
    cerr << "The executable should be invoked with exactly one filepath, C, M, "
            "Q1, Q2, Q3 "
            "argument. Example ./MyImageApplication '../../Lena_512_512.rgb' 1 "
            "1 8 8 8"
         << endl;
    exit(1);
  }
  cout << "First argument: " << wxApp::argv[0] << endl;
  cout << "Second argument: " << wxApp::argv[1] << endl;
  cout << "C: " << wxApp::argv[2] << endl;
  cout << "M: " << wxApp::argv[3] << endl;
  cout << "Q1: " << wxApp::argv[4] << endl;
  cout << "Q2: " << wxApp::argv[5] << endl;
  cout << "Q3: " << wxApp::argv[6] << endl;

  string imagePath = wxApp::argv[1].ToStdString();
  int C = stoi(wxApp::argv[2].ToStdString());
  int M = stoi(wxApp::argv[3].ToStdString());
  int Q1 = stoi(wxApp::argv[4].ToStdString());
  int Q2 = stoi(wxApp::argv[5].ToStdString());
  int Q3 = stoi(wxApp::argv[6].ToStdString());

  if (Q1 < 1 || Q1 > 8 || Q2 < 1 || Q2 > 8 || Q3 < 1 || Q3 > 8) {
    cerr << "Q1, Q2, Q3 must be between 1 and 8" << endl;
    exit(1);
  }

  if (C != 1 && C != 2 && C != 0) {
    cerr << "C must be either 2 or 1" << endl;
    exit(1);
  }

  if (M != 1 && M != 2) {
    cerr << "M must be either 2 or 1" << endl;
    exit(1);
  }
  MyFrame *frame = new MyFrame("Image Display", imagePath, C, M, Q1, Q2, Q3);
  frame->Show(true);

  // return true to continue, false to exit the application
  return true;
}

/**
 * Constructor for the MyFrame class.
 * Here we read the pixel data from the file and set up the scrollable window.
 */
MyFrame::MyFrame(const wxString &title, string imagePath, int C, int M, int Q1,
                 int Q2, int Q3)
    : wxFrame(NULL, wxID_ANY, title) {

  // Modify the height and width values here to read and display an image with
  // different dimensions.
  width = 352;
  height = 288;
  unsigned char *outData;
  unsigned char *inData = readImageData(imagePath, width, height);
  if (C == 0) {
    testCases();
    exit(0);
  }

  if (C == 1 | C == 2) {
    outData = quantize(inData, width, height, C, M, Q1, Q2, Q3);

  } else {
    outData = inData;
  }

  // the last argument is static_data, if it is false, after this call the
  // pointer to the data is owned by the wxImage object, which will be
  // responsible for deleting it. So this means that you should not delete the
  // data yourself.
  inImage.SetData(inData, width, height, false);

  outImage.SetData(outData, width, height, false);

  // Set up the scrolled window as a child of this frame
  scrolledWindow = new wxScrolledWindow(this, wxID_ANY);
  scrolledWindow->SetScrollbars(10, 10, width * 2, height);
  scrolledWindow->SetVirtualSize(width * 2, height);

  // Bind the paint event to the OnPaint function of the scrolled window
  scrolledWindow->Bind(wxEVT_PAINT, &MyFrame::OnPaint, this);

  // Set the frame size
  SetClientSize(width * 2, height);

  // Set the frame background color
  SetBackgroundColour(*wxBLACK);
}

/**
 * The OnPaint handler that paints the UI.
 * Here we paint the image pixels into the scrollable window.
 */
void MyFrame::OnPaint(wxPaintEvent &event) {
  wxBufferedPaintDC dc(scrolledWindow);
  scrolledWindow->DoPrepareDC(dc);

  wxBitmap inImageBitmap = wxBitmap(inImage);
  dc.DrawBitmap(inImageBitmap, 0, 0, false);

  wxBitmap outImageBitmap = wxBitmap(outImage);
  dc.DrawBitmap(outImageBitmap, width, 0, false);
}

double *matrixMultiply(double *A, int hA, int wA, double *B, int hB, int wB) {
  if (wA != hB) {
    throw invalid_argument("Incompatible matrix dimensions for multiplication");
  }
  

  double *C = new double[hA * wB];


  for (int i = 0; i < hA; i++) {
    for (int j = 0; j < wB; j++) {
      double sum = 0;
      for (int k = 0; k < wA; k++) {
        double a = A[i * wA + k];
        double b = B[k * wB + j];

        sum += a * b;
      }
      C[i * wB + j] = static_cast<double>(sum);

    }
  }
#if DEBUG

#endif
  return C;
}

double *from255to1(unsigned char *imageData, int width, int height) {
  double *outData = (double *)malloc(width * height * 3 * sizeof(double));
  if (!outData)
    return NULL;
#if DEBUG
  cout << "Converting from [0,255] to [0,1]" << endl;
#endif
  for (int i = 0; i < height * width; i++) {
    outData[3 * i] = (double)imageData[3 * i] / 255.0;
    outData[3 * i + 1] = (double)imageData[3 * i + 1] / 255.0;
    outData[3 * i + 2] = (double)imageData[3 * i + 2] / 255.0;
  }
#if DEBUG
  for (int i = 0; i < 10; i++) {
    cout << "[from255to1] Pixel " << i << ": " << (int)imageData[3 * i]
         << " -> " << (double)outData[3 * i] << endl;
  }
#endif

  return outData;
}

unsigned char *from1to255(double *imageData, int width, int height) {
  unsigned char *outData =
      (unsigned char *)malloc(width * height * 3 * sizeof(unsigned char));
  if (!outData)
    return NULL;
#if DEBUG
  cout << "Converting from [0,1] to [0,255]" << endl;
#endif
  for (int i = 0; i < height * width; i++) {
    if (imageData[3 * i] > 2) {
#if DEBUG
      cout << "Error: pixel value greater than 1: " << (int)imageData[3 * i]
           << endl;
#endif
      return NULL;
    }
    // clamp to [0,1] before scaling to prevent wrap-around artifacts
    double r = imageData[3 * i];
    double g = imageData[3 * i + 1];
    double b = imageData[3 * i + 2];
    if (r < 0.0) r = 0.0; else if (r > 1.0) r = 1.0;
    if (g < 0.0) g = 0.0; else if (g > 1.0) g = 1.0;
    if (b < 0.0) b = 0.0; else if (b > 1.0) b = 1.0;
    outData[3 * i] = (int)(r * 255.0);
    outData[3 * i + 1] = (int)(g * 255.0);
    outData[3 * i + 2] = (int)(b * 255.0);
  }
#if DEBUG
  for (int i = 0; i < 10; i++) {
    cout << "[from1to255] Pixel " << i << ": " << (double)imageData[3 * i]
         << " -> " << (int)outData[3 * i] << endl;
  }
#endif
  return outData;
}

unsigned char *rgbrgbrgb2rrggbb(unsigned char *imageData, int width,
                                int height) {
  unsigned char *outData =
      (unsigned char *)malloc(width * height * 3 * sizeof(unsigned char));
  if (!outData)
    return NULL;
  for (int i = 0; i < height * width; i++) {
    outData[3 * i] = imageData[3 * i];
    outData[3 * i + 1] = imageData[3 * i + 1];
    outData[3 * i + 2] = imageData[3 * i + 2];
  }
  return outData;
}

unsigned char *rrggbb2rgbrgbrgb(unsigned char *imageData, int width,
                                int height) {
  unsigned char *outData =
      (unsigned char *)malloc(width * height * 3 * sizeof(unsigned char));
  if (!outData)
    return NULL;
  for (int i = 0; i < height * width; i++) {
    outData[3 * i] = imageData[3 * i];
    outData[3 * i + 1] = imageData[3 * i + 1];
    outData[3 * i + 2] = imageData[3 * i + 2];
  }
  return outData;
}

double *RGB2YUV(unsigned char *imageData, int width, int height) {
  unsigned char *outData =
      (unsigned char *)malloc(width * height * 3 * sizeof(unsigned char));
  if (!outData)
    return NULL;
#if DEBUG
  cout << "Converting from RGB to YUV" << endl;
#endif
  double *outDataDouble = (double *)malloc(width * height * 3 * sizeof(double));
  if (!outDataDouble)
    return NULL;

  double *imageDataDouble =
      (double *)malloc(width * height * 3 * sizeof(double));
  if (!imageDataDouble)
    return NULL;

  // BT.601 RGB to YUV conversion matrix
  double RGB2YUVMatrix[] = {0.299, 0.587, 0.114,  -0.147, -0.289,
                            0.436, 0.615, -0.515, -0.100};

  imageDataDouble = from255to1(imageData, width, height);
#if DEBUG
  double minY=1e9, maxY=-1e9, minU=1e9, maxU=-1e9, minV=1e9, maxV=-1e9;
  for (int i=0;i<width*height;i++) {
      minY = std::min(minY, imageDataDouble[3*i]);
      maxY = std::max(maxY, imageDataDouble[3*i]);
      minU = std::min(minU, imageDataDouble[3*i+1]);
      maxU = std::max(maxU, imageDataDouble[3*i+1]);
      minV = std::min(minV, imageDataDouble[3*i+2]);
      maxV = std::max(maxV, imageDataDouble[3*i+2]);
  }
  cout << "Y:[" << minY << "," << maxY << "] "
      << "U:[" << minU << "," << maxU << "] "
      << "V:[" << minV << "," << maxV << "]" << endl;

#endif
  for (int i = 0; i < height * width; i++) {
    // We populate YUV values of each pixel in that order
    // YUV.YUV.YUV and so on for all pixels
    double r = imageDataDouble[3 * i];
    double g = imageDataDouble[3 * i + 1];
    double b = imageDataDouble[3 * i + 2];
    
    // Ensure RGB values are in valid range [0,1]
    if (r < 0.0) r = 0.0; else if (r > 1.0) r = 1.0;
    if (g < 0.0) g = 0.0; else if (g > 1.0) g = 1.0;
    if (b < 0.0) b = 0.0; else if (b > 1.0) b = 1.0;

    double B[] = {r,g,b};

    double *yuv =
        matrixMultiply(RGB2YUVMatrix, 3, 3, B, 3, 1);
#if DEBUG
    if (i < 10) {
      cout << "[RGB2YUV] Pixel " << i << ": R=" << r << " G=" << g << " B=" << b
           << " => Y=" << yuv[0] << " U=" << yuv[1] << " V=" << yuv[2] << endl;
    }
#endif
    outDataDouble[3 * i] = yuv[0];
    outDataDouble[3 * i + 1] = yuv[1];
    outDataDouble[3 * i + 2] = yuv[2];
    delete[] yuv;
  }
#if DEBUG
  minY=1e9, maxY=-1e9, minU=1e9, maxU=-1e9, minV=1e9, maxV=-1e9;
  for (int i=0;i<width*height;i++) {
      minY = std::min(minY, outDataDouble[3*i]);
      maxY = std::max(maxY, outDataDouble[3*i]);
      minU = std::min(minU, outDataDouble[3*i+1]);
      maxU = std::max(maxU, outDataDouble[3*i+1]);
      minV = std::min(minV, outDataDouble[3*i+2]);
      maxV = std::max(maxV, outDataDouble[3*i+2]);
  }
  cout << "Y:[" << minY << "," << maxY << "] "
      << "U:[" << minU << "," << maxU << "] "
      << "V:[" << minV << "," << maxV << "]" << endl;

#endif
  // outData = from1to256(outDataDouble, width, height);
  // delete [] imageDataDouble;
  // delete [] outDataDouble;
  return outDataDouble;
}

unsigned char *YUV2RGB(double *imageData, int width, int height) {
  unsigned char *outData =
      (unsigned char *)malloc(width * height * 3 * sizeof(unsigned char));
  if (!outData)
    return NULL;

  double *outDataDouble = (double *)malloc(width * height * 3 * sizeof(double));
  if (!outDataDouble)
    return NULL;

// double *imageDataDouble =
//     (double *)malloc(width * height * 3 * sizeof(double));
// if (!imageDataDouble) return NULL;
#if DEBUG
  cout << "Converting from YUV to RGB" << endl;
#endif
  // imageDataDouble = from256to1(imageData, width, height);
  // BT.601 YUV to RGB conversion matrix (inverse of RGB2YUV)
  double YUV2RGBMatrix[] = {1.0,     0.0, 1.1398, 1.0, -0.3946,
                            -0.5806, 1.0, 2.0321, 0.0};
  for (int i = 0; i < height * width; i++) {
    // We populate RGB values of each pixel in that order
    // RGB.RGB.RGB and so on for all pixels

    double y = imageData[3 * i];
    double u = imageData[3 * i + 1];
    double v = imageData[3 * i + 2];
    
    // Clamp YUV values to valid ranges to prevent artifacts
    if (y < 0.0) y = 0.0; else if (y > 1.0) y = 1.0;
    if (u < -0.436) u = -0.436; else if (u > 0.436) u = 0.436;
    if (v < -0.615) v = -0.615; else if (v > 0.615) v = 0.615;
    
#if DEBUG
    if (i < 10) {
      cout << "[YUV2RGB] Clamped Y=" << y << " U=" << u << " V=" << v << endl;
    }
#endif
    double B[] = {y,u,v};

    double *rgb =
        matrixMultiply(YUV2RGBMatrix, 3, 3, B, 3, 1);
#if DEBUG
    if (i < 10) {
      cout << "[YUV2RGB] Pixel " << i << ": Y=" << y << " U=" << u << " V=" << v
           << " => R=" << rgb[0] << " G=" << rgb[1] << " B=" << rgb[2] << endl;
    }
#endif

    // clamp RGB to [0,1] to avoid invalid values after inverse transform
    double r = rgb[0];
    double g = rgb[1];
    double b = rgb[2];
    if (r < 0.0) r = 0.0; else if (r > 1.0) r = 1.0;
    if (g < 0.0) g = 0.0; else if (g > 1.0) g = 1.0;
    if (b < 0.0) b = 0.0; else if (b > 1.0) b = 1.0;
    outDataDouble[3 * i] = r;
    outDataDouble[3 * i + 1] = g;
    outDataDouble[3 * i + 2] = b;
    delete[] rgb;
  }
  outData = from1to255(outDataDouble, width, height);
  // delete [] imageDataDouble;
  // delete [] outDataDouble;

  return outData;
}

double *nonUniformQuantizeYUV(const double *imageData, int width, int height,
                              int Q1, int Q2, int Q3, double Ymin, double Ymax,
                              double Umin, double Umax, double Vmin,
                              double Vmax, int H) {
#if DEBUG
  cout << "[nonUniformQuantizeYUV] Non-uniform quantization with C=2 Q1=" << Q1 << " Q2=" << Q2 << " Q3=" << Q3 << endl;
#endif

  const int size = width * height;
  if (!imageData || size <= 0 || H <= 1)
    return nullptr;

  auto clampi = [](int x, int lo, int hi) {
    return x < lo ? lo : (x > hi ? hi : x);
  };

  const double lo[3] = {Ymin, Umin, Vmin};
  const double hi[3] = {Ymax, Umax, Vmax};
  const double span[3] = {Ymax - Ymin, Umax - Umin, Vmax - Vmin};

  
  if (span[0] <= 0 || span[1] <= 0 || span[2] <= 0)
    return nullptr;


  const int levels[3] = {1 << Q1, 1 << Q2, 1 << Q3};
  if (levels[0] <= 1 || levels[1] <= 1 || levels[2] <= 1)
    return nullptr;
  double *out = (double *)std::malloc(size * 3 * sizeof(double));
  if (!out)
    return nullptr;

  for (int c = 0; c < 3; ++c) {
    std::vector<int> hist(H, 0);
    for (int i = 0; i < size; ++i) {
      double v = imageData[3 * i + c];
      if (!(v == v))
        v = lo[c];
      if (v < lo[c])
        v = lo[c];
      if (v > hi[c])
        v = hi[c];

      int bin = (int)std::floor((v - lo[c]) / span[c] * (H - 1));
      bin = clampi(bin, 0, H - 1);
      hist[bin]++;
    }
#if DEBUG
    cout << "[nonUniformQuantizeYUV] Channel " << c << " histogram: ";
    for (int i = 0; i < H; ++i) {
      if (hist[i] > 0)
        cout << "[" << i << "]=" << hist[i] << " ";
    }
    cout << endl;
#endif
    // Lloyd-Max
    const int L = levels[c];
    std::vector<int> boundaryBin(L + 1, 0);
    std::vector<double> rep(L, 0.0);

    {
      std::vector<int> cdf(H, 0);
      cdf[0] = hist[0];
      for (int i = 1; i < H; ++i)
        cdf[i] = cdf[i - 1] + hist[i];
      boundaryBin[0] = 0;
      int curLevel = 1;
      for (int i = 0; i < H && curLevel < L; ++i) {
        if ((long long)cdf[i] * L >= (long long)curLevel * size) {
          boundaryBin[curLevel++] = i;
        }
      }
      boundaryBin[L] = H - 1;
    }

    for (int k = 0; k < L; ++k) {
      double a = lo[c] + (boundaryBin[k] + 0.5) * (span[c] / H);
      double b = lo[c] + (boundaryBin[k + 1] + 0.5) * (span[c] / H);
      rep[k] = 0.5 * (a + b);
    }

    const int maxIter = 25;
    const double eps = 1e-6;
    for (int iter = 0; iter < maxIter; ++iter) {
      double maxDelta = 0.0;
      for (int k = 0; k < L; ++k) {
        long long wsum = 0;
        long long cnt = 0;
        int b0 = boundaryBin[k];
        int b1 = boundaryBin[k + 1];
        for (int b = b0; b <= b1; ++b) {
          long long h = hist[b];
          if (!h)
            continue;
          double v = lo[c] + (b + 0.5) * (span[c] / H);
          wsum += (long long)(h);
          cnt += 1;
        }
        double newRep;
        if (wsum == 0) {
          newRep = rep[k];
        } else {
          long double num = 0.0L;
          long double den = 0.0L;
          for (int b = b0; b <= b1; ++b) {
            if (hist[b] == 0)
              continue;
            long double h = (long double)hist[b];
            long double v = (long double)(lo[c] + (b + 0.5) * (span[c] / H));
            num += h * v;
            den += h;
          }
          newRep = (double)(num / (den > 0 ? den : 1.0L));
        }
        maxDelta = std::max(maxDelta, std::fabs(newRep - rep[k]));
        rep[k] = newRep;
      }

      for (int k = 1; k < L; ++k) {
        double mid = 0.5 * (rep[k - 1] + rep[k]);
        int binIdx = (int)std::floor((mid - lo[c]) / span[c] * H - 0.5);
        binIdx = clampi(binIdx, 0, H - 1);
        if (binIdx <= boundaryBin[k - 1])
          binIdx = boundaryBin[k - 1] + 1;
        if (binIdx >= boundaryBin[k + 1])
          binIdx = boundaryBin[k + 1] - 1;
        boundaryBin[k] = binIdx;
      }

      if (maxDelta < eps)
        break;
    }

    std::vector<int> bin2lvl(H, 0);
    int lvl = 0;
    for (int i = 0; i < H; ++i) {
      while (lvl + 1 < L && i > boundaryBin[lvl + 1])
        ++lvl;
      bin2lvl[i] = lvl;
    }

    for (int i = 0; i < size; ++i) {
      double v = imageData[3 * i + c];
      if (!(v == v))
        v = lo[c];
      if (v < lo[c])
        v = lo[c];
      if (v > hi[c])
        v = hi[c];

      int bin = (int)std::floor((v - lo[c]) / span[c] * (H - 1));
      bin = clampi(bin, 0, H - 1);
      out[3 * i + c] = rep[bin2lvl[bin]];
    }
  }

  return out;
}

unsigned char *nonUniformQuantize(unsigned char *imageData, int width,
                                  int height, int C, int Q1, int Q2, int Q3) {

  unsigned char *outData = (unsigned char *)malloc(width * height * 3);
  if (!outData)
    return NULL;

#if DEBUG
  cout << "[nonUniformQuantize] Non-uniform quantization with C=" << C
       << " Q1=" << Q1 << " Q2=" << Q2 << " Q3=" << Q3 << endl;
#endif
  int Rlevels = 1 << Q1;
  int Glevels = 1 << Q2;
  int Blevels = 1 << Q3;
  int levels[3] = {Rlevels, Glevels, Blevels};
  int size = width * height;

  for (int c = 0; c < 3; c++) { // 0=R, 1=G, 2=B
    int hist[256] = {0};

    for (int i = 0; i < size; i++) {
      hist[imageData[3 * i + c]]++;
    }

    int cdf[256] = {0};
    cdf[0] = hist[0];
    for (int i = 1; i < 256; i++)
      cdf[i] = cdf[i - 1] + hist[i];

    int boundaries[levels[c] + 1];
    boundaries[0] = 0;
    int curLevel = 1;
    for (int i = 0; i < 256 && curLevel < levels[c]; i++) {
      if (cdf[i] >= curLevel * size / levels[c]) {
        boundaries[curLevel++] = i;
      }
    }
    boundaries[levels[c]] = 255;

#if DEBUG
    cout << "Boundaries for channel " << c << ": ";
    for (int i = 0; i <= levels[c]; i++) {
      cout << boundaries[i] << " ";
    }
    cout << endl;
#endif

    unsigned char rep[levels[c]];
    for (int i = 0; i < levels[c]; i++) {
      rep[i] = (boundaries[i] + boundaries[i + 1]) / 2;
    }

    for (int i = 0; i < size; i++) {
      unsigned char val = imageData[3 * i + c];
      for (int l = 0; l < levels[c]; l++) {
        if (val >= boundaries[l] && val <= boundaries[l + 1]) {
          outData[3 * i + c] = rep[l];
          break;
        }
      }
    }
  }

  return outData;
}

unsigned char *quantize(unsigned char *imageData, int width, int height, int C,
                        int M, int Q1, int Q2, int Q3) {
  unsigned char *outData = (unsigned char *)malloc(width * height * 3);
  if (!outData)
    return NULL;

  double *imageDataDouble =
      (double *)malloc(width * height * 3 * sizeof(double));
  if (!imageDataDouble)
    return NULL;

  double *outDataDouble = (double *)malloc(width * height * 3 * sizeof(double));
  if (!outDataDouble)
    return NULL;

  if (C == 2) {
    imageDataDouble = RGB2YUV(imageData, width, height);
    if (M == 1) {
      const int YRange = 1 << Q1; // N_Y
      const int URange = 1 << Q2; // N_U
      const int VRange = 1 << Q3; // N_V

      // YUV ranges based on BT.601 standard (should match actual conversion output)
      const double Ya = 0.0, Yb = 1.0;
      const double Ua = -0.436, Ub = 0.436;
      const double Va = -0.615, Vb = 0.615;

      const double Y_span = Yb - Ya;
      const double U_span = Ub - Ua;
      const double V_span = Vb - Va;

      for (int i = 0; i < width * height; ++i) {
        double y = imageDataDouble[3 * i];
        double u = imageDataDouble[3 * i + 1];
        double v = imageDataDouble[3 * i + 2];

        // ---- Y ----
        int yLevel = (int)floor((y - Ya) / Y_span * YRange);
        if (yLevel < 0)
          yLevel = 0;
        else if (yLevel >= YRange)
          yLevel = YRange - 1;
        double yInterval = Y_span / YRange;
        double yq = Ya + (yLevel + 0.5) * yInterval;

        // ---- U ----
        int uLevel = (int)floor((u - Ua) / U_span * URange);
        if (uLevel < 0)
          uLevel = 0;
        else if (uLevel >= URange)
          uLevel = URange - 1;
        double uInterval = U_span / URange;
        double uq = Ua + (uLevel + 0.5) * uInterval;

        // ---- V ----
        int vLevel = (int)floor((v - Va) / V_span * VRange);
        if (vLevel < 0)
          vLevel = 0;
        else if (vLevel >= VRange)
          vLevel = VRange - 1;
        double vInterval = V_span / VRange;
        double vq = Va + (vLevel + 0.5) * vInterval;

        outDataDouble[3 * i] = yq;
        outDataDouble[3 * i + 1] = uq;
        outDataDouble[3 * i + 2] = vq;
      }
    } else if (M == 2) {
#if DEBUG
      cout << "imageDataDouble size: " << width * height * 3 << endl;
#endif

      outDataDouble =
          nonUniformQuantizeYUV(imageDataDouble, width, height, Q1, Q2, Q3);
      if (!outDataDouble) {
        cout << "Error in non-uniform quantization" << endl;
        free(outData);
        return NULL;
      }


    } else {
      free(outData);
      return NULL;
    }
#if DEBUG
    cout << "Probe for outDataDouble after quantization:" << endl;
    for (int i = 0; i < 10; i++) {
      cout << "[quantize] Pixel " << i << ": Y=" << outDataDouble[3 * i]
           << " U=" << outDataDouble[3 * i + 1]
           << " V=" << outDataDouble[3 * i + 2] << endl;
    }
#endif
    outData = YUV2RGB(outDataDouble, width, height);

  } else {
    if (M == 1) {
      unsigned int RRange = 1 << Q1;
      unsigned int GRange = 1 << Q2;
      unsigned int BRange = 1 << Q3;

      unsigned int rInterval = 256 / RRange;
      unsigned int gInterval = 256 / GRange;
      unsigned int bInterval = 256 / BRange;

      for (int i = 0; i < width * height; i++) {
        unsigned int rLevel = (imageData[3 * i] * RRange) / 256;
        unsigned int gLevel = (imageData[3 * i + 1] * GRange) / 256;
        unsigned int bLevel = (imageData[3 * i + 2] * BRange) / 256;

        if (rLevel >= RRange)
          rLevel = RRange - 1;
        if (gLevel >= GRange)
          gLevel = GRange - 1;
        if (bLevel >= BRange)
          bLevel = BRange - 1;

        outData[3 * i] = rLevel * rInterval + rInterval / 2;
        outData[3 * i + 1] = gLevel * gInterval + gInterval / 2;
        outData[3 * i + 2] = bLevel * bInterval + bInterval / 2;
      }

    } else if (M == 2) {
      outData = nonUniformQuantize(imageData, width, height, C, Q1, Q2, Q3);
    } else {
      free(outData);
      return NULL;
    }
  }

  return outData;
}

unsigned int getPictureError(unsigned char *img1, unsigned char *img2, int width,
                        int height) {
  if (!img1 || !img2 || width <= 0 || height <= 0)
    return -1;

  unsigned int error = 0;

  for (int i = 0; i < width * height * 3; i++) {
    int diff = pow(abs((int)img1[i] - (int)img2[i]),2);
    error += diff;
  }
  error = error / (width * height * 3);
  return error;
}

struct TestResult {
  int N;
  int C, M, Q1, Q2, Q3;
  unsigned int error;
  
  TestResult(int n, int c, int m, int q1, int q2, int q3, unsigned int err) 
    : N(n), C(c), M(m), Q1(q1), Q2(q2), Q3(q3), error(err) {}
};

void testCases() {
  string imagePath = "./mountain_352x288.rgb";
  int width = 352;
  int height = 288;
  unsigned char *inData = readImageData(imagePath, width, height);
  unsigned char *outData;

  vector<TestResult> results;
  unsigned int minError = UINT_MAX;
  string minMode;
  int N[3] = {4,6,8};
  
  for (int C=1; C<=2; C++) {
    for (int M=1; M<=2; M++) {
      for (int i=0; i<3; i++) {
        for (int j=1; j<=N[i]; j++) {
          for (int k=N[i]-j-1; k>0; k--) {
            if (1){
              int l = N[i] - j - k;
              int Q1 = j;
              int Q2 = k;
              int Q3 = l;
              outData = quantize(inData, width, height, C, M, Q1, Q2, Q3);
              unsigned int error = getPictureError(inData, outData, width, height);
              
              // Store result
              results.push_back(TestResult(Q1+Q2+Q3, C, M, Q1, Q2, Q3, error));
              
              if (error < minError) {
                minError = error;
                minMode = "C=" + to_string(C) + " M=" + to_string(M) + " Q1=" + to_string(Q1) + " Q2=" + to_string(Q2) + " Q3=" + to_string(Q3);
              }
              cout << "Testing C=" << C << " M=" << M << " Q1=" << Q1 << " Q2=" << Q2 << " Q3=" << Q3 << " N="<< Q1+Q2+Q3 << endl;
              cout << "Total absolute error for <"<< C << ", "<< M<<", "<<Q1<<", "<<Q2<<", "<<Q3<<">: "<< error << endl;

              string mode = (C==1) ? "RGB" : "YUV";
              mode += (M==1) ? "_equal_" : "_nonuniform_";
              string fileNameBase = "out_" + mode +"N"+to_string(N[i])+"_"+ to_string(Q1) + to_string(Q2) + to_string(Q3)+"_error_"+to_string(error);
              // saveToRgbFile(fileNameBase + ".rgb", outData, width, height);
              // saveToPngFile(fileNameBase + ".png", outData, width, height);
            }
          }
        }
      }
    }
  }
  
  // Sort results by N in ascending order
  sort(results.begin(), results.end(), [](const TestResult& a, const TestResult& b) {
    return a.N < b.N;
  });
  
  // Save results to table file
  ofstream tableFile("error_table.txt");
  if (tableFile.is_open()) {
    tableFile << "N\t<C,M,Q1,Q2,Q3>\t\tError" << endl;
    for (const auto& result : results) {
      tableFile << result.N << "\t<" << result.C << "," << result.M << "," 
                << result.Q1 << "," << result.Q2 << "," << result.Q3 << ">\t\t" 
                << result.error << endl;
    }
    tableFile.close();
    cout << "Error table saved to error_table.txt" << endl;
  } else {
    cerr << "Error: Could not create error_table.txt" << endl;
  }
  
  cout << "Minimum error found: " << minError << endl;
  cout << "Achieved with mode: " << minMode << endl;

  // outData = quantize(inData, width, height, 1, 1, 8, 8, 8);
  // saveToRgbFile("out_RGB_equal_888.rgb", outData, width, height);
  // saveToPngFile("out_RGB_equal_888.png", outData, width, height);

  // outData = quantize(inData, width, height, 1, 1, 2, 2, 2);
  // saveToRgbFile("out_RGB_equal_222.rgb", outData, width, height);
  // saveToPngFile("out_RGB_equal_222.png", outData, width, height);

  // outData = quantize(inData, width, height, 2, 2, 2, 3, 3);
  // saveToRgbFile("out_YUV_opt_233.rgb", outData, width, height);
  // saveToPngFile("out_YUV_opt_233.png", outData, width, height);
}

void saveToPngFile(string filePath, unsigned char *data, int width,
                   int height) {
  wxImage image;
  image.SetData(data, width, height, false);
  if (!image.SaveFile(filePath, wxBITMAP_TYPE_PNG)) {
    cerr << "Error Saving PNG File" << endl;
    exit(1);
  }
}

void saveToRgbFile(string filePath, unsigned char *data, int width,
                   int height) {
  ofstream outputFile(filePath, ios::binary);
  if (!outputFile.is_open()) {
    cerr << "Error Opening File for Writing" << endl;
    exit(1);
  }
  vector<char> Rbuf(width * height);
  vector<char> Gbuf(width * height);
  vector<char> Bbuf(width * height);

  for (int i = 0; i < width * height; i++) {
    Rbuf[i] = data[3 * i];
    Gbuf[i] = data[3 * i + 1];
    Bbuf[i] = data[3 * i + 2];
  }

  outputFile.write(Rbuf.data(), width * height);
  outputFile.write(Gbuf.data(), width * height);
  outputFile.write(Bbuf.data(), width * height);

  outputFile.close();
}

/** Utility function to read image data */
unsigned char *readImageData(string imagePath, int width, int height) {

  // Open the file in binary mode
  ifstream inputFile(imagePath, ios::binary);

  if (!inputFile.is_open()) {
    cerr << "Error Opening File for Reading" << endl;
    exit(1);
  }

  // Create and populate RGB buffers
  vector<char> Rbuf(width * height);
  vector<char> Gbuf(width * height);
  vector<char> Bbuf(width * height);

  /**
   * The input RGB file is formatted as RRRR.....GGGG....BBBB.
   * i.e the R values of all the pixels followed by the G values
   * of all the pixels followed by the B values of all pixels.
   * Hence we read the data in that order.
   */

  inputFile.read(Rbuf.data(), width * height);
  inputFile.read(Gbuf.data(), width * height);
  inputFile.read(Bbuf.data(), width * height);

  inputFile.close();

  /**
   * Allocate a buffer to store the pixel values
   * The data must be allocated with malloc(), NOT with operator new. wxWidgets
   * library requires this.
   */
  unsigned char *inData =
      (unsigned char *)malloc(width * height * 3 * sizeof(unsigned char));

  for (int i = 0; i < height * width; i++) {
    // We populate RGB values of each pixel in that order
    // RGB.RGB.RGB and so on for all pixels
    inData[3 * i] = Rbuf[i];
    inData[3 * i + 1] = Gbuf[i];
    inData[3 * i + 2] = Bbuf[i];
  }

  return inData;
}

wxIMPLEMENT_APP(MyApp);