#include <cmath>
#include <wx/wx.h>
#include <wx/dcbuffer.h>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#define DEBUG 1
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
  MyFrame(const wxString &title, string imagePath, int C, int M, int Q1, int Q2, int Q3);

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
double *from256to1(unsigned char * imageData, int width, int height);
unsigned char *from1to256(double * imageData, int width, int height);
unsigned char *quantize(unsigned char * imageData, int width, int height,int C, int M, int Q1, int Q2, int Q3);
double *RGB2YUV(unsigned char * imageData, int width, int height);
unsigned char *YUV2RGB(double * imageData, int width, int height);
void saveToPngFile(string filePath, unsigned char *data, int width, int height);
void saveToRgbFile(string filePath, unsigned char *data, int width, int height);
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
    cerr << "The executable should be invoked with exactly one filepath, C, M, Q1, Q2, Q3 "
            "argument. Example ./MyImageApplication '../../Lena_512_512.rgb' 1 1 8 8 8"
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
MyFrame::MyFrame(const wxString &title, string imagePath, int C, int M, int Q1, int Q2, int Q3)
    : wxFrame(NULL, wxID_ANY, title) {

  // Modify the height and width values here to read and display an image with
  // different dimensions.    
  width = 512;
  height = 512;
  unsigned char *outData;
  unsigned char *inData = readImageData(imagePath, width, height);
  if (C==0){
    testCases();
    exit(0);
  }
      
  if(C == 1 | C == 2){
    outData = quantize(inData, width, height, C, M, Q1, Q2, Q3);

  }else{
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
  scrolledWindow->SetScrollbars(10, 10, width*2, height);
  scrolledWindow->SetVirtualSize(width*2, height);

  // Bind the paint event to the OnPaint function of the scrolled window
  scrolledWindow->Bind(wxEVT_PAINT, &MyFrame::OnPaint, this);

  // Set the frame size
  SetClientSize(width*2, height);

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

double* matrixMultiply(
    double* A, int hA, int wA,
    double* B, int hB, int wB
) {
    if (wA != hB) {
        throw invalid_argument("Incompatible matrix dimensions for multiplication");
    }

    double* C = new double[hA * wB];

    #if DEBUG
      // cout << "Parameters: hA:" << hA << " wA:" << wA << " hB:" << hB << " wB:" << wB << endl; 
      // for (int i = 0; i < hA; i++) {
          // for (int j = 0; j < wA; j++) {
              // cout << A[i * wA + j] << " ";
          // }
          // cout << endl;
      // }
    #endif
    for (int i = 0; i < hA; i++) {
        for (int j = 0; j < wB; j++) {
            double sum = 0;
            for (int k = 0; k < wA; k++) {
                double a = A[i * wA + k];
                double b = B[k * wB + j];
                #if DEBUG
                  // if (i == 0 && j < 10) {
                  //     cout << "A[0][" << k << "] = " << a << ", B[" << k << "][" << j << "] = " << b << endl;
                  // }
                #endif
                sum += a * b;
            }
            #if DEBUG
              // if(sum > 255) {
              //     cout << "pixel value overflow: " << sum << endl;
              // }else if (sum < 0){
              //     cout << "pixel value underflow: " << sum << endl;
              // }
            #endif
            // if (sum > 255) sum = 255;
            // if (sum < 0) sum = 0;
            C[i * wB + j] = static_cast<double>(sum);
            #if DEBUG
              // if (i == 0 && j < 10) {
              //     cout << "C[0][" << j << "] = " << C[i * wB + j] << endl;
              // }
            #endif
        }
    }
    return C;
}

double *from256to1(unsigned char * imageData, int width, int height){
    double *outData =
      (double *)malloc(width * height * 3 * sizeof(double));
    if (!outData) return NULL;
    #if DEBUG
      cout << "Converting from [0,255] to [0,1]" << endl;
    #endif
    for (int i = 0; i < height * width; i++) {
        outData[3*i] = (double)imageData[3*i] / 255.0;
        outData[3*i+1] = (double)imageData[3*i] / 255.0;
        outData[3*i+2] = (double)imageData[3*i] / 255.0;
    }
    #if DEBUG
      for (int i = 0; i < 10; i++) {
          cout << "Pixel " << i << ": " << (int)imageData[3*i] << " -> " << (double)outData[3*i] << endl;
      }
    #endif
    
    return outData;
}

unsigned char *from1to256(double * imageData, int width, int height){
    unsigned char *outData =
      (unsigned char *)malloc(width * height * 3 * sizeof(unsigned char));
    if (!outData) return NULL;
    #if DEBUG
      cout << "Converting from [0,1] to [0,255]" << endl;
    #endif
    for (int i = 0; i < height * width; i++) {
        if(imageData[3*i] > 2){
          #if DEBUG
            cout << "Error: pixel value greater than 1: " << (int)imageData[3*i] << endl;
          #endif
          return NULL;
        } 
        outData[3*i] = (int)(imageData[3*i] * 255);
        outData[3*i+1] = (int)(imageData[3*i] * 255);
        outData[3*i+2] = (int)(imageData[3*i] * 255);
    }
    return outData;
}

unsigned char *rgbrgbrgb2rrggbb(unsigned char * imageData, int width, int height){
    unsigned char *outData =
      (unsigned char *)malloc(width * height * 3 * sizeof(unsigned char));
    if (!outData) return NULL;
    for (int i = 0; i < height * width; i++) {
        outData[3*i] = imageData[3*i];
        outData[3*i+1] = imageData[3*i];
        outData[3*i+2] = imageData[3*i];
    }
    return outData;
}

unsigned char *rrggbb2rgbrgbrgb(unsigned char * imageData, int width, int height){
    unsigned char *outData =
      (unsigned char *)malloc(width * height * 3 * sizeof(unsigned char));
    if (!outData) return NULL;
    for (int i = 0; i < height * width; i++) {
        outData[3*i] = imageData[3*i];
        outData[3*i+1] = imageData[3*i+1];
        outData[3*i+2] = imageData[3*i+2];
    }
    return outData;
}

double *RGB2YUV(unsigned char * imageData, int width, int height) {
  unsigned char *outData =
      (unsigned char *)malloc(width * height * 3 * sizeof(unsigned char));
  if (!outData) return NULL;
  #if DEBUG
    cout << "Converting from RGB to YUV" << endl;
  #endif
  double *outDataDouble =
      (double *)malloc(width * height * 3 * sizeof(double));
  if (!outDataDouble) return NULL;
  
  double *imageDataDouble =
      (double *)malloc(width * height * 3 * sizeof(double));
  if (!imageDataDouble) return NULL;

  

  double RGB2YUVMatrix[] = {
      0.299, 0.587, 0.114,
      -0.147, -0.289, 0.436,
      0.615, -0.515, -0.100
  };

  imageDataDouble = from256to1(imageData, width, height);

  for (int i = 0; i < height * width; i++) {
    // We populate YUV values of each pixel in that order
    // YUV.YUV.YUV and so on for all pixels
    double r = imageDataDouble[3*i];
    double g = imageDataDouble[3*i + 1];
    double b = imageDataDouble[3*i + 2];
    double *yuv = matrixMultiply(RGB2YUVMatrix, 3, 3, (double[]){r, g, b}, 3, 1);
    #if DEBUG
      if (i < 10) {
          cout << "Pixel " << i << ": R=" << r << " G=" << g << " B=" << b 
               << " => Y=" << yuv[0] << " U=" << yuv[1] << " V=" << yuv[2] << endl;
      }
    #endif
    outDataDouble[3 * i] = yuv[0];
    outDataDouble[3 * i + 1] = yuv[1];
    outDataDouble[3 * i + 2] = yuv[2];
    delete[] yuv;
    
  }
  // outData = from1to256(outDataDouble, width, height);
  // delete [] imageDataDouble;
  // delete [] outDataDouble;
  return outDataDouble;
  
}

unsigned char *YUV2RGB(double * imageData, int width, int height) {
  unsigned char *outData =
      (unsigned char *)malloc(width * height * 3 * sizeof(unsigned char));
  if (!outData) return NULL;

  double *outDataDouble =
      (double *)malloc(width * height * 3 * sizeof(double));
  if (!outDataDouble) return NULL;
  
  // double *imageDataDouble =
  //     (double *)malloc(width * height * 3 * sizeof(double));
  // if (!imageDataDouble) return NULL;
  #if DEBUG
    cout << "Converting from YUV to RGB" << endl;
  #endif
  // imageDataDouble = from256to1(imageData, width, height);
  double YUV2RGBMatrix[] = {
      1.0, 0.0, 1.1398,
      1.0, -0.3946, -0.5806,
      1.0, 2.0321, 0.0
  };
  for (int i = 0; i < height * width; i++) {
    // We populate RGB values of each pixel in that order
    // RGB.RGB.RGB and so on for all pixels
    double y = imageData[3*i];
    double u = imageData[3*i + 1];
    double v = imageData[3*i + 2];
    double *rgb = matrixMultiply(YUV2RGBMatrix, 3, 3, (double[]){y, u, v}, 3, 1);
    outDataDouble[3 * i] = rgb[0];
    outDataDouble[3 * i + 1] = rgb[1];
    outDataDouble[3 * i + 2] = rgb[2];
    delete[] rgb;
    
  }
  outData = from1to256(outDataDouble, width, height);
  // delete [] imageDataDouble;
  // delete [] outDataDouble;

  
  return outData;
  
}
unsigned char *nonUniformQuantize(unsigned char* imageData, int width, int height, int Q1, int Q2, int Q3){
    int Rlevels = 1 << Q1;      
    int Glevels = 1 << Q2;
    int Blevels = 1 << Q3;     
    int levels[3] = {Rlevels, Glevels, Blevels};
    int size = width * height;
    unsigned char *outData = (unsigned char*)malloc(width * height * 3);
    if (!outData) return NULL;

    for (int c = 0; c < 3; c++) { // 0=R, 1=G, 2=B
        int hist[256] = {0};

        for (int i = 0; i < size; i++) {
            hist[imageData[3*i + c]]++;
        }


        int cdf[256] = {0};
        cdf[0] = hist[0];
        for (int i = 1; i < 256; i++) cdf[i] = cdf[i-1] + hist[i];


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
            rep[i] = (boundaries[i] + boundaries[i+1]) / 2;
        }


        for (int i = 0; i < size; i++) {
            unsigned char val = imageData[3*i + c];
            for (int l = 0; l < levels[c]; l++) {
                if (val >= boundaries[l] && val <= boundaries[l+1]) {
                    outData[3*i + c] = rep[l];
                    break;
                }
            }
        }
    }

    return outData;
}


unsigned char *quantize(unsigned char * imageData, int width, int height, int C, int M, int Q1, int Q2, int Q3) {
  unsigned char *outData = (unsigned char *)malloc(width * height * 3);
    if (!outData) return NULL;
    
    double *imageDataDouble =
      (double *)malloc(width * height * 3 * sizeof(double));
    if (!imageDataDouble) return NULL;

    double *outDataDouble =
      (double *)malloc(width * height * 3 * sizeof(double));
    if (!outDataDouble) return NULL;

    if (C == 2){
      imageDataDouble = RGB2YUV(imageData, width, height);
      if (M == 1){
        const int YRange = 1 << Q1;  // N_Y
        const int URange = 1 << Q2;  // N_U
        const int VRange = 1 << Q3;  // N_V

        const double Ya = 0.0,     Yb = 1.0;      
        const double Ua = -0.436,  Ub = 0.436;   
        const double Va = -0.615,  Vb = 0.615; 

        const double Y_span = Yb - Ya;
        const double U_span = Ub - Ua;
        const double V_span = Vb - Va;

        for (int i = 0; i < width * height; ++i) {
            double y = imageDataDouble[3*i];
            double u = imageDataDouble[3*i+1];
            double v = imageDataDouble[3*i+2];

            // ---- Y ----
            int yLevel = (int)floor((y - Ya) / Y_span * YRange);
            if (yLevel < 0) yLevel = 0;
            else if (yLevel >= YRange) yLevel = YRange - 1;
            double yInterval = Y_span / YRange;
            double yq = Ya + (yLevel + 0.5) * yInterval;

            // ---- U ----
            int uLevel = (int)floor((u - Ua) / U_span * URange);
            if (uLevel < 0) uLevel = 0;
            else if (uLevel >= URange) uLevel = URange - 1;
            double uInterval = U_span / URange;
            double uq = Ua + (uLevel + 0.5) * uInterval;

            // ---- V ----
            int vLevel = (int)floor((v - Va) / V_span * VRange);
            if (vLevel < 0) vLevel = 0;
            else if (vLevel >= VRange) vLevel = VRange - 1;
            double vInterval = V_span / VRange;
            double vq = Va + (vLevel + 0.5) * vInterval;

            outDataDouble[3*i]   = yq;
            outDataDouble[3*i+1] = uq;
            outDataDouble[3*i+2] = vq;
        }  
      } else if (M == 2){
        outData = nonUniformQuantize(imageData, width, height, Q1, Q2, Q3);
      }else{
        free(outData);
        return NULL;
      }
      outData = YUV2RGB(outDataDouble, width, height);

    }else{
      if (M == 1){
      unsigned int RRange = 1 << Q1;
      unsigned int GRange = 1 << Q2;
      unsigned int BRange = 1 << Q3;

      unsigned int rInterval = 256 / RRange;
      unsigned int gInterval = 256 / GRange;
      unsigned int bInterval = 256 / BRange;

      for (int i = 0; i < width * height; i++) {
          unsigned int rLevel = (imageData[3*i] * RRange) / 256;
          unsigned int gLevel = (imageData[3*i+1] * GRange) / 256;
          unsigned int bLevel = (imageData[3*i+2] * BRange) / 256;

          if (rLevel >= RRange) rLevel = RRange - 1;
          if (gLevel >= GRange) gLevel = GRange - 1;
          if (bLevel >= BRange) bLevel = BRange - 1;

          outData[3*i] = rLevel * rInterval + rInterval / 2;
          outData[3*i+1] = gLevel * gInterval + gInterval / 2;
          outData[3*i+2] = bLevel * bInterval + bInterval / 2;
      }

      
      

      } else if (M == 2){
        outData = nonUniformQuantize(imageData, width, height, Q1, Q2, Q3);
        
        
      }else{
        free(outData);
        return NULL;
      }
    }

    return outData;
  
}



void testCases(){
  string imagePath = "../../Lena_512_512.rgb";
  int width = 512;
  int height = 512;
  unsigned char *inData = readImageData(imagePath, width, height);
  unsigned char *outData;

  outData = quantize(inData, width, height,1, 1, 8, 8, 8);
  saveToRgbFile("out_RGB_equal_888.rgb", outData, width, height);
  saveToPngFile("out_RGB_equal_888.png", outData, width, height);

  outData = quantize(inData, width, height,1, 1, 2, 2, 2);
  saveToRgbFile("out_RGB_equal_222.rgb", outData, width, height);
  saveToPngFile("out_RGB_equal_222.png", outData, width, height);

  outData = quantize(inData, width, height,2, 2, 2, 3, 3);
  saveToRgbFile("out_YUV_opt_233.rgb", outData, width, height);
  saveToPngFile("out_YUV_opt_233.png", outData, width, height);
  

}

void saveToPngFile(string filePath, unsigned char *data, int width, int height){
    wxImage image;
    image.SetData(data, width, height, false);
    if (!image.SaveFile(filePath, wxBITMAP_TYPE_PNG)) {
        cerr << "Error Saving PNG File" << endl;
        exit(1);
    }

}

void saveToRgbFile(string filePath, unsigned char *data, int width, int height){
    ofstream outputFile(filePath, ios::binary);
    if (!outputFile.is_open()) {
        cerr << "Error Opening File for Writing" << endl;
        exit(1);
    }
    vector<char> Rbuf(width * height);
    vector<char> Gbuf(width * height);
    vector<char> Bbuf(width * height);

    for (int i = 0; i < width * height; i++) {
        Rbuf[i] = data[3*i];
        Gbuf[i] = data[3*i + 1];
        Bbuf[i] = data[3*i + 2];
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