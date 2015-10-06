#include "stdafx.h"
#include <omp.h>
#include "ImgFacotryAPI.h"


ImgFacotryAPI::ImgFacotryAPI(unsigned int width, unsigned int height, unsigned int singlebyte, string path) :_width(width), _height(height), _singlebyte(singlebyte), _path(path)
{
    cv::initModule_nonfree();
    _imgcontent = nullptr;

    _imgcontent = new unsigned int[width * height];
    SetCoordinates();
}

void ImgFacotryAPI::LonLat2MapXY(const double& lon, const double& lat, const double& height, double& mapx_, double& mapy_)
{
    double lon_[1], lat_[1], height_[1];
    lon_[0] = lon; lat_[0] = lat; height_[0] = height;
    poTransform->Transform(1, lon_, lat_, height_);
    mapx_ = lon_[0];
    mapy_ = lat_[0];
}


bool SelfBuild(MatchToolType& _lefttool, VirtualMachine* _leftModel)
{
    if (_lefttool._typeofFac == Enum_Match_RigurousModel)
    {
        _lefttool._ptool = (_leftModel);
    }
    else if (_lefttool._typeofFac == Enum_Match_GEOTIFF)
    {
        GDALDataset* poDataset = (GDALDataset*)GDALOpen(_lefttool.filepath.c_str(), GA_ReadOnly);
        _lefttool._ptool = new double[6];
        double transtemp[6];
        int xedge = poDataset->GetRasterXSize() - 1;
        int yedge = poDataset->GetRasterYSize() - 1;
        poDataset->GetGeoTransform(transtemp);
        ///计算反过来的投影信息
        double realtrans[2][3];
        Mat trans_mat(2, 3, CV_64F, &realtrans);
        Point2f dst_pf[3] = { Point2f(0, 0), Point2f(xedge, 0), Point2f(0, yedge) };
        Point2f src_pf[3] = { Point2f(transtemp[0], transtemp[3]),
            Point2f(transtemp[0] + transtemp[1] * xedge, transtemp[3] + transtemp[4] * xedge),
            Point2f(transtemp[0] + transtemp[2] * yedge, transtemp[3] + transtemp[5] * yedge) };

        trans_mat = getAffineTransform(src_pf, dst_pf);

        for (int i = 0; i < 2; i++)
        for (int j = 0; j < 3; j++)
        {
            static_cast<double*>(_lefttool._ptool)[i * 3 + j] = trans_mat.at<double>(i, j);
        }
        GDALClose(poDataset);
    }
    else if (_lefttool._typeofFac == Enum_Match_IMGHEADER)
    {
        GDALDataset* poDataset = (GDALDataset*)GDALOpen(_lefttool.filepath.c_str(), GA_ReadOnly);
        _lefttool._ptool = new double[6];
        poDataset->GetGeoTransform(static_cast<double*>(_lefttool._ptool));
        GDALClose(poDataset);
    }

    return true;
}

bool PreMatch(MatchToolType& _lefttool, MatchToolType& _right, VirtualMachine* _leftModel, int nonsense)
{
    return true;
}

void ImgFacotryAPI::Translate(const MatchToolType& _lefttool, const MatchToolType& _rightTool, const double& height, double* transx, double* transy, int sizeall)
{
    double loninfi, latinfo;
    for (int i = 0; i < sizeall; i++)
    {
        _lefttool.FromXY2LONLAT(transx[i], transy[i], height, loninfi, latinfo);
        double mapx, mapy;
        LonLat2MapXY(loninfi, latinfo, height, mapx, mapy);
        _rightTool.FromRect2XY(mapx, mapy, height, transx[i], transy[i]);
    }
}

using namespace std;
void ImgFacotryAPI::Begin_Match_WorkFlow(string srcImg, string refimg, VirtualMachine* info)
{
    lsm_multi_times = true;
    MatchToolType _left;
    _left.filepath = srcImg;
    _left.outputpath = srcImg;
    _left._typeofFac = Enum_Match_RigurousModel;

    MatchToolType _right;
    _right.filepath = refimg;
    _right.outputpath = _right.filepath + "PreMatch.tiff";
    _right._typeofFac = Enum_Match_GEOTIFF;

    cout << "跳过预处理?" << endl;
    int mode = -1;
    cin >> mode;
    if (mode == 0)
        PreMatch(_left, _right, info, nullptr);

    int pausefor_histo = 0;
    cout << "进行直方图匹配:??" << endl;
    cin >> pausefor_histo;
    string oldpath = _left.outputpath;
    _left.outputpath = "D:\\OpenCVWithCuda\\RPCModel\\ConsoleApplication1\\ConsoleApplication1\\Whole.tif21_.tif";

    if (pausefor_histo)
    {
        Histo_Match_Image(_left, _right);
    }

    _left.outputpath = oldpath;
    string oldright = _right.outputpath;
    _right.outputpath = _right.outputpath + "8bit.tif";


    cout << "直方图均衡化完毕 ，请将8位后的均衡化影像放到路径" << _right.outputpath.c_str() << endl;
    int pausefor_preparation;
    cin >> pausefor_preparation;
    Match_Pixel(_left, _right);

    _left.outputpath = "D:\\OpenCVWithCuda\\RPCModel\\ConsoleApplication1\\ConsoleApplication1\\Whole.tif25_.tif";
    _right.outputpath = oldright;// "D:\\IRS\\JB11-1_IRS_000116063_003_002_L0\\JB11-1_IRS_000116063_003_004_002_L3C_TOABT.tiffPreMatch.tiff";
    vector<DMatch> emptyMatch;

    cout << "跳过相位匹配" << endl;
    int phase_mode = 0;
    cin >> phase_mode;
    if (!phase_mode)
        Match_Phase_Correlation(_left, _right, emptyMatch);

    LSM_Match(_left, _right, 19);
    cout << "" << endl;
}

#include <opencv2/nonfree/nonfree.hpp>
///粗匹配结果上进行小于0.1个像素精度的匹配
void ImgFacotryAPI::Match_Phase_Correlation(MatchToolType& _left, MatchToolType& _right, vector<DMatch> matchres)
{
    //GDALDataset* poData_set_Left = static_cast<GDALDataset*>(GDALOpen(_left.outputpath.c_str(), GA_ReadOnly));
    //GDALDataset* poData_set_Right = static_cast<GDALDataset*>(GDALOpen(_right.outputpath.c_str(), GA_ReadOnly));
    //GDALRasterBand* poBand_Left = poData_set_Left->GetRasterBand(1);
    //GDALRasterBand* poBand_Right = poData_set_Right->GetRasterBand(1);
    Mat poData_set_Left_Ori = imread(_left.outputpath, CV_LOAD_IMAGE_ANYDEPTH);
    Mat poData_set_Right = imread(_right.outputpath, CV_LOAD_IMAGE_ANYDEPTH);
    Mat poData_set_Left; poData_set_Left_Ori.convertTo(poData_set_Left, CV_64F);

    int left_width = poData_set_Left.cols;
    int left_hei = poData_set_Left.rows;

    int rigt_width = poData_set_Right.cols;
    int rigt_hei = poData_set_Right.rows;

    float* left_image_data = new float[left_width * left_hei];
    double* rigt_image_data = new double[rigt_width * rigt_hei];
    memset(left_image_data, 0, sizeof(float)* left_width * left_hei);
    memset(rigt_image_data, 0, sizeof(double)* rigt_width * rigt_hei);



    ///获取匹配点
    int matchsize = _leftkey_.size();

    ///初始搜索窗口大小为16 * 16个像素
    int windowsize = 33;
    Mat left_Window_Info;
    Mat right_Window_Info;
    int Loc_Min_Left_X = 0;
    int Loc_Min_Left_Y = 0;
    int Loc_Min_Right_X = 0;
    int Loc_Min_Right_Y = 0;
    ///对每一个匹配点过滤
    heivaluecount = 0;

    ofstream fftwofs;
    fftwofs.open("C:\\Temp\\fftwshift.xls");
    for (int i = 0; i < matchsize; i++)
    {
        int half_windows_size = windowsize / 2;
        Loc_Min_Left_X = _leftkey_[i].x - half_windows_size;
        Loc_Min_Left_Y = _leftkey_[i].y - half_windows_size;
        Loc_Min_Right_X = _rightkey_[i].x - half_windows_size;
        Loc_Min_Right_Y = _rightkey_[i].y - half_windows_size;
        if (Loc_Min_Left_X < 0 || Loc_Min_Left_Y < 0 || Loc_Min_Right_X < 0 || Loc_Min_Right_Y < 0)
        {
            cout << "Out of Range! Skipped" << endl;
            continue;
        }
        int Loc_Max_Left_X = _leftkey_[i].x + half_windows_size;
        int Loc_Max_Left_Y = _leftkey_[i].y + half_windows_size;
        int Loc_Max_Right_X = _rightkey_[i].x + half_windows_size;
        int Loc_Max_Right_Y = _rightkey_[i].y + half_windows_size;
        if (Loc_Max_Left_X > left_width - 1 || Loc_Max_Left_Y > left_hei - 1 || Loc_Max_Right_X > rigt_width - 1 || Loc_Max_Right_Y > rigt_hei - 1)
        {
            cout << "Out of Range! Skipped" << endl;
            continue;
        }

        left_Window_Info = poData_set_Left(Rect(Loc_Min_Left_X, Loc_Min_Left_Y, windowsize, windowsize));
        right_Window_Info = poData_set_Right(Rect(Loc_Min_Right_X, Loc_Min_Right_Y, windowsize, windowsize));

        double maxvalue = 0;
        Point2d pointmove;
        is_shifted = false;
        Phase_Correction(left_Window_Info, right_Window_Info, maxvalue, pointmove);
        fftwofs.precision(30);
        if (true)
        {
            double x_left = _leftkey_[i].x;
            double y_left = _leftkey_[i].y;
            fftwofs << x_left << " " << 10785 - y_left << " " << _rightkey_[i].x << " " << 10785 - _rightkey_[i].y << endl;
        }
        ///避免边缘点
        //cout << maxvalue;

    }

    //GDALClose(poData_set_Left);
    //GDALClose(poData_set_Right);
}

void printParams(cv::Algorithm* algorithm) {
    std::vector<std::string> parameters;
    algorithm->getParams(parameters);

    for (int i = 0; i < (int)parameters.size(); i++) {
        std::string param = parameters[i];
        int type = algorithm->paramType(param);
        std::string helpText = algorithm->paramHelp(param);
        std::string typeText;

        switch (type) {
        case cv::Param::BOOLEAN:
            typeText = "bool";
            break;
        case cv::Param::INT:
            typeText = "int";
            break;
        case cv::Param::REAL:
            typeText = "real (double)";
            break;
        case cv::Param::STRING:
            typeText = "string";
            break;
        case cv::Param::MAT:
            typeText = "Mat";
            break;
        case cv::Param::ALGORITHM:
            typeText = "Algorithm";
            break;
        case cv::Param::MAT_VECTOR:
            typeText = "Mat vector";
            break;
        }
        std::cout << "Parameter '" << param << "' type=" << typeText << " help=" << helpText << std::endl;
    }
}

void ImgFacotryAPI::Match_Pixel(MatchToolType& _lefttool, MatchToolType& _rightool)
{
    initModule_nonfree();


    string leftpath = _lefttool.outputpath;
    string rightpath = _rightool.outputpath;

    Mat img_left = imread(leftpath.c_str(), CV_LOAD_IMAGE_GRAYSCALE);
    Mat img_right = imread(rightpath.c_str(), CV_LOAD_IMAGE_GRAYSCALE);


    cv::Ptr<cv::FeatureDetector> detector = cv::FeatureDetector::create("PyramidFAST");
    //printParams(detector);

    //detector->set("edgeThreshold", 1);
    cv::Ptr<cv::DescriptorExtractor> descriptor = cv::DescriptorExtractor::create("ORB");

    //detector->
    ///获得关键点
    std::vector<cv::KeyPoint> keypoints1, keypoints2, keypionts3;


    //imshow("hei",picu8bit);
    //waitKey(0);
    cv::initModule_nonfree();
    //detector->detect(picu8bit(rectleft), keypionts3);


    cv::initModule_nonfree();
    //detector->detect(img_1, keypoints1);
    detector->detect(img_left, keypoints1);
    detector->detect(img_right, keypoints2);
    cout << "关键点数目： k1:" << keypoints1.size() << "," << keypoints2.size() << endl;
    ///extract features
    Mat descriptions_1, descriptions_2;

    descriptor->compute(img_left, keypoints1, descriptions_1);
    descriptor->compute(img_right, keypoints2, descriptions_2);

    BFMatcher matcher(NORM_HAMMING, true);
    vector<DMatch> matches1;
    matcher.match(descriptions_1, descriptions_2, matches1);

    double max_dist = 0; double min_dist = 9999;

    ///Calculation
    double sumall = 0;
    double avg_angle = 0;
    double avg_all = 0;
    int counter = 0;
    vector<double> anglevect(matches1.size());
    for (int i = 0; i < matches1.size(); i++)
    {
        double dist = matches1[i].distance;
        double index_x_1 = keypoints1[matches1[i].queryIdx].pt.x;
        double index_y_1 = keypoints1[matches1[i].queryIdx].pt.y;
        double index_x_2 = keypoints2[matches1[i].trainIdx].pt.x;
        double index_y_2 = keypoints2[matches1[i].trainIdx].pt.y;
        double directionangle = atan(fabs((index_y_2 - index_y_1) / (fabs(index_x_2 - index_x_1) + 480))) * 180 / CV_PI;
        if (directionangle < 20)
        {
            counter++;
            avg_angle += fabs(directionangle);
        }
        anglevect[i] = fabs(directionangle);
        min_dist = (dist < min_dist) ? dist : min_dist;
        max_dist = (dist > max_dist) ? dist : max_dist;
        avg_all += fabs(dist) / (matches1.size() - 1);

    }

    avg_angle /= counter;
    double rms_error_angle = 0;
    for (int i = 0; i < matches1.size(); i++)
    {
        double dist = matches1[i].distance;

        min_dist = (dist < min_dist) ? dist : min_dist;
        max_dist = (dist > max_dist) ? dist : max_dist;
        if (anglevect[i] < 20)
            rms_error_angle += pow(anglevect[i] - avg_angle, 2.0) / (counter - 1);

        sumall += pow(dist - avg_all, 2.0) / (matches1.size() - 1);
    }
    sumall = sqrt(sumall);
    rms_error_angle = sqrt(rms_error_angle);
    cout.precision(20);
    cout << "Match size is:" << matches1.size() << endl;
    cout << "avg dst is " << avg_all << endl;
    cout << "Rms error is" << sumall << endl;
    cout << "----Mat dist : " << max_dist << endl;
    cout << "----Min dist : " << min_dist << endl;
    cout << "Avg angle is" << avg_angle << endl;
    cout << "Angle rms angle is" << rms_error_angle << endl;
    std::vector<DMatch> good_matches;

    ofstream ofs;
    for (int i = 0; i < matches1.size(); i++)
    {
        if (fabs(anglevect[i] - avg_angle) > 2 * rms_error_angle)
            continue;


        if (fabs(matches1[i].distance - avg_all) < sumall * 3)
        {
            double index_x_1 = keypoints1[matches1[i].queryIdx].pt.x;
            double index_y_1 = keypoints1[matches1[i].queryIdx].pt.y;
            double index_x_2 = keypoints2[matches1[i].trainIdx].pt.x;
            double index_y_2 = keypoints2[matches1[i].trainIdx].pt.y;
            Point2d hei_left(index_x_1, index_y_1);
            Point2d hei_right(index_x_2, index_y_2);
            _leftkey_.push_back(hei_left);
            _rightkey_.push_back(hei_right);
            good_matches.push_back(matches1[i]);
        }

    }

    cout << "Pruned match: " << good_matches.size() << endl;
    Mat img_matches;
    drawMatches(img_left, keypoints1, img_right, keypoints2,
        good_matches, img_matches, Scalar::all(-1), Scalar::all(-1),
        vector<char>(), DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);
    imwrite("d:\\matchresult.tif", img_matches);
}

double ImgFacotryAPI::CalcCorrelation(double * pl, double *pr, int wm, int wn)
{
    double sum_r = 0, sum_l = 0;			//左右灰度和
    double sum_r_2 = 0, sum_l_2 = 0;		//左右灰度平方和
    double sum_lr = 0;					//左右灰度成绩和
    int i, j;		//循环变量
    for (j = 0; j < wn; j++)
    {
        for (i = 0; i < wm; i++)
        {
            int index_id = i + 1 + (j + 1) * (wm + 2);
            sum_l += pl[index_id];
            sum_r += pr[index_id];
            sum_l_2 += pl[index_id] * pl[index_id];
            sum_r_2 += pr[index_id] * pr[index_id];
            sum_lr += pl[index_id] * pr[index_id];
        }
    }
    if ((sum_l_2 - sum_l*sum_l / wm / wn)*(sum_r_2 - sum_r*sum_r / wm / wn) == 0)
    {
        return 1;
    }
    double re = (sum_lr - sum_l*sum_r / wm / wn) / sqrt((sum_l_2 - sum_l*sum_l / wm / wn)*(sum_r_2 - sum_r*sum_r / wm / wn));
    return re;

}

bool ImgFacotryAPI::Histo_Match_Image(MatchToolType& leftinfi, MatchToolType& rightinfi)
{

    Mat leftMat_Ori = imread(leftinfi.outputpath, CV_LOAD_IMAGE_ANYDEPTH);
    Mat leftMat; leftMat_Ori.convertTo(leftMat, CV_64F);
    Mat rightMat_Ori = imread(leftinfi.outputpath, CV_LOAD_IMAGE_ANYDEPTH);

    int bins = static_cast<int>(pow(2, 12));
    int hist_size[] = { bins };
    float range[] = { 0, static_cast<float>(bins) };
    const float* ranges[] = { range };
    Mat heiinfo;

    MatND hiswhole; int channels[] = { 0 };
    /* calcHist( &gray, 1, channels, Mat(), // do not use mask
    hist, 1, hist_size, ranges,
    true, // the histogram is uniform
    false );  */
    //int hei = info->at<unsigned short>(1000,5);
    Mat leftMat_Int; leftMat_Ori.convertTo(leftMat_Int, CV_16U);
    Mat rightMat_Int;  rightMat_Ori.convertTo(rightMat_Int, CV_16U);
    calcHist(&leftMat_Int, 1, channels, Mat(), hiswhole, 1, hist_size, ranges, true, false);
    //DrawHist(hiswhole);

    
    
    float* look_uptable = new float[4096];
    HistoMatching(rightinfi.outputpath, hiswhole, 12, look_uptable);

    //后面的就交给GDAL
    GDALDriver* poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    int leftwidth = leftMat_Ori.cols;
    int leftheight = leftMat_Ori.rows;
    GDALDataset* poDatasetOri = (GDALDataset*)GDALOpen(rightinfi.outputpath.c_str(), GA_ReadOnly);
    GDALRasterBand* poBand_Ori = poDatasetOri->GetRasterBand(1);
    GDALDataset* poDataset = poDriver->Create((rightinfi.outputpath +"hei.tif" ).c_str(), leftwidth, leftheight, 1, GDT_Float64, NULL);
    double* dataoutput = new double[leftwidth * leftheight];
    GDALRasterBand* poBand = poDataset->GetRasterBand(1);
    poBand_Ori->RasterIO(GF_Read, 0, 0, leftwidth, leftheight, dataoutput, leftwidth, leftheight, GDT_Float64, 0, 0);
    for (int i = 0; i < leftwidth * leftheight; i++)
    {
        int heivalue = static_cast<int>(dataoutput[i]);
        dataoutput[i] = look_uptable[heivalue];
    }

    poBand->RasterIO(GF_Write, 0, 0, leftwidth, leftheight, dataoutput, leftwidth, leftheight, GDT_Float64, 0, 0);
    GDALClose(poDataset);
    GDALClose(poDatasetOri);
    return true;
}


bool ImgFacotryAPI::LSM_Match(MatchToolType leftinfi, MatchToolType rightinfi, const vector<Point2d>& leftkey_, const vector<Point2d>& rightkey_,
    vector<Point2d>& _leftnewkey, vector<Point2d>& _rightnewkey,
    vector<MatchResult>& _lsmresult,int windowsize)
{
    CSatOrbit helper;
    if (windowsize % 2 == 0)
    {
        cout << "请将最小二乘窗口大小设置为奇数" << endl;
        return false;
    }

    int halfwindow = (windowsize) / 2;
    Mat leftMat_Ori = imread(leftinfi.outputpath, CV_LOAD_IMAGE_ANYDEPTH);
    Mat leftMat; leftMat_Ori.convertTo(leftMat, CV_64F);
    Mat rightMat_Ori = imread(rightinfi.outputpath, CV_LOAD_IMAGE_ANYDEPTH);


    Mat rightMat = rightMat_Ori;
    int leftwidth_ = leftMat.cols;
    int leftheight_ = leftMat.rows;
    int rightwidth_ = rightMat.cols;
    int rightheight_ = rightMat.rows;  ///checked

    Rect leftrange(0, 0, leftwidth_, leftheight_);
    Rect rightrange(0, 0, rightwidth_, rightheight_);
    vector<Point2d> _left_lsm;
    vector<Point2d> _right_lsm;
    vector<Point2d> _hfactor;
    vector<Point2d> _leftmove;
    vector<Point2d> _rightmove;
    vector<double> _correvaluearray;
    ///对每个匹配点进行最小二乘精匹配
    int icc = 0;
    int ncount_count_valid = 0;
    for (int ii = 0; ii < leftkey_.size(); ii++)
    {
        cout << ii << "\\" << leftkey_.size() << "\\" << ncount_count_valid << "\r";
        ResampleKernels* my_kernel = new ResampleKernels;
        my_kernel->Init(Raised_Cosine_6P);
        int realwidth_ = windowsize - 2;
        int realheight_ = windowsize - 2;
        Rect arealeft(leftkey_[ii].x - halfwindow + 1, leftkey_[ii].y - halfwindow + 1, windowsize, windowsize);
        Rect arearight(rightkey_[ii].x - halfwindow + 1, rightkey_[ii].y - halfwindow + 1, windowsize, windowsize);
        if (!safeopencv::Inside(leftrange, arealeft)) continue; /// 此点在边界，暂时被忽略

        ///系数矩阵
        double* coeffarray = new double[realwidth_ * realwidth_ * 8];
        double* staticarray = new double[realwidth_ * realheight_]; // 观测值
        double* tarray = new double[realwidth_ * realheight_ * 8];  // 常量矩阵

        double CTC[64] = { 0 };
        double CTL[8] = { 0 };

        double pt[4] = { 0 };
        double P0[4] = { leftkey_[ii].x, leftkey_[ii].y, rightkey_[ii].x, rightkey_[ii].y };

        ///左右片变形参数
        double a0 = 0, a1 = 1, a2 = 0, b0 = 0, b1 = 0, b2 = 1;
        ///辐射畸变参数
        double h0 = 0;
        double h1 = 1;

        ///左影像数据、指针
        Mat Left_Block_ori = leftMat(arealeft).clone();
        Mat Left_Block; Left_Block_ori.convertTo(Left_Block, CV_64F);
        void* pvoidleft = static_cast<void*>(Left_Block.data);
        double* pimgleft = new double[windowsize * windowsize];// = static_cast<double*>(pvoidleft);

        ///左影像像平面坐标,真特么麻烦
        vector<Point2f> Left_Block_Pixel;
        safeopencv::Mat_Fill(Left_Block_Pixel, windowsize, Enum_Mat_Fill_Symmetric);

        ///填充左影像数据
        for (int i = 0; i < windowsize * windowsize; i++)
        {
            Point2f leftresample_pos(Left_Block_Pixel[i].x + leftkey_[ii].x, Left_Block_Pixel[i].y + leftkey_[ii].y);
            double tempvalue = safeopencv::Resample(leftMat, leftresample_pos, my_kernel);
            pimgleft[i] = tempvalue;
        }

        ///左、右窗口。  左窗口不动,右窗口初始值为粗匹配计算结果
        ///开始最小二乘参数的解算
        double re[2] = { 0.95, -1 };
        double th0, th1, ta0, ta1, ta2, tb1, tb2, tb0;
        int ncounter = 0;
        do
        {
            if (ncounter++ > 20)
                break;
            //判断右窗口是否仍在影像范围内
            if (!safeopencv::Inside(rightrange, arearight)) {
                cout << "右窗口溢出！！" << endl;
            }
            vector<Point2f> Right_Blcok_Pixel(windowsize * windowsize);
            double a_affine_trans_data[6] = { a1, a2, a0,
                b1, b2, b0 };
            Mat affine_mat(2, 3, CV_32F, &a_affine_trans_data);
            Point2f pixelstart;
            pixelstart.x = rightkey_[ii].x; pixelstart.y = rightkey_[ii].y;
            safeopencv::Pixel_Warp_Affine(Left_Block_Pixel, Right_Blcok_Pixel, pixelstart, a_affine_trans_data);


            ///重采样, 辐射畸变改正
            Mat Right_Block(windowsize, windowsize, CV_64F);
            int window = windowsize * windowsize;
            void* pvoidright = static_cast<void*>(Right_Block.data);
            double* pimgright = static_cast<double*>(pvoidright);

            for (int q = 0; q < window; q++)
            {
                //cout << safeopencv::Resample(rightMat, Right_Blcok_Pixel[q], my_kernel) << endl;
                double tempvalue = safeopencv::Resample(rightMat, Right_Blcok_Pixel[q], my_kernel);
                pimgright[q] = h0 + h1 * tempvalue;
            }
            ///Right_Block_Pixel 为了方便重采样,此时是原影响坐标,需要转换
            ///计算相关系数
            re[1] = CalcCorrelation(pimgleft, pimgright, windowsize - 2, windowsize - 2);
            //cout << re[1] << endl;

            if (re[1]<re[0] && re[0]>0.95)
                break;

            re[0] = re[1];

            double a[8], aa[64], al[8];
            double x[8];
            cout.precision(4);
            Mat hei_test_left = leftMat(arealeft);
            Mat hei_test_right = rightMat(arearight);
            for (int j = 0; j < windowsize - 2; j++)
            {
                for (int i = 0; i < windowsize - 2; i++)
                {
                    int index_now = (j * (windowsize - 2) + i) * 8;
                    coeffarray[index_now + 0] = 1;
                    coeffarray[index_now + 1] = pimgleft[(j + 1) * windowsize + i + 1];
                    coeffarray[index_now + 2] = h1*(pimgright[(j + 1) * windowsize + i + 2] - pimgright[(j + 1) * windowsize + i])*0.5;
                    coeffarray[index_now + 3] = h1*(pimgright[(j + 1) * windowsize + i + 2] - pimgright[(j + 1) * windowsize + i])*0.5*Left_Block_Pixel[(j + 1) * windowsize + i + 1].x;//lx[i + 1][j + 1];
                    coeffarray[index_now + 4] = h1*(pimgright[(j + 1) * windowsize + i + 2] - pimgright[(j + 1) * windowsize + i])*0.5*Left_Block_Pixel[(j + 1) * windowsize + i + 1].y;//ly[i + 1][j + 1];
                    coeffarray[index_now + 5] = h1*(pimgright[(j + 2) * windowsize + i + 1] - pimgright[(j + 0) * windowsize + i + 1])*0.5;
                    coeffarray[index_now + 6] = h1*(pimgright[(j + 2) * windowsize + i + 1] - pimgright[(j + 0) * windowsize + i + 1])*0.5*Left_Block_Pixel[(j + 1) * windowsize + i + 1].x;
                    coeffarray[index_now + 7] = h1*(pimgright[(j + 2) * windowsize + i + 1] - pimgright[(j + 0) * windowsize + i + 1])*0.5*Left_Block_Pixel[(j + 1) * windowsize + i + 1].y;

                    staticarray[j*(windowsize - 2) + i] = pimgleft[(j + 1) * windowsize + i + 1] - pimgright[(j + 1) * windowsize + i + 1];



                }

            }


            helper.transpose(coeffarray, tarray, realwidth_ * realheight_, 8);
            helper.mult(tarray, coeffarray, CTC, 8, realwidth_*realheight_, 8);
            helper.mult(tarray, staticarray, CTL, 8, realwidth_*realheight_, 1);



            //CTL为8X1
            if (!helper.invers_matrix(CTC, 8))		 					//求逆
            {
                cout << "不可逆" << endl;
                ;
            }

            double V[8];


            helper.mult(CTC, CTL, V, 8, 8, 1);
            //memcpy(x, al, sizeof(double)* 8);
            cout.precision(20);


            //V对应dh0 dh1  da0-da2 db0-db2
            th0 = h0; th1 = h1; ta0 = a0; ta1 = a1; ta2 = a2; tb1 = b1; tb2 = b2; tb0 = b0;
            h0 = th0 + V[0] + th0*V[1];
            h1 = th1 + th1*V[1];
            a0 = ta0 + V[2] + ta0*V[3] + tb0*V[4];
            a1 = ta1 + ta1*V[3] + tb1*V[4];
            a2 = ta2 + ta2*V[3] + tb2*V[4];
            b0 = tb0 + V[5] + ta0*V[6] + tb0*V[7];
            b1 = tb1 + ta1*V[6] + tb1*V[7];
            b2 = tb2 + ta2*V[6] + tb2*V[7];


        } while (true);


        double sum_x = 0, sum_y = 0, sum_gx = 0, sum_gy = 0;
        for (int j = 0; j < realheight_; j++)
        {
            for (int i = 0; i < realwidth_; i++)
            {
                sum_x += Left_Block_Pixel[i + 1 + (j + 1) *windowsize].x * (pimgleft[i + 2 + (j + 1) * windowsize] - pimgleft[i + (j + 1) * windowsize])*(pimgleft[i + 2 + (j + 1) * windowsize] - pimgleft[i + (j + 1) * windowsize])*0.25;
                sum_y += Left_Block_Pixel[i + 1 + (j + 1) *windowsize].y * (pimgleft[i + 1 + (j + 2) * windowsize] - pimgleft[i + 1 + j * windowsize])*(pimgleft[i + 1 + (j + 2) * windowsize] - pimgleft[i + 1 + j * windowsize])*0.25;
                sum_gx += (pimgleft[i + 2 + (j + 1) * windowsize] - pimgleft[i + (j + 1) * windowsize])*(pimgleft[i + 2 + (j + 1) * windowsize] - pimgleft[i + (j + 1) * windowsize])*0.25;
                sum_gy += (pimgleft[i + 1 + (j + 2) * windowsize] - pimgleft[i + 1 + j * windowsize])*(pimgleft[i + 1 + (j + 2) * windowsize] - pimgleft[i + 1 + j * windowsize])*0.25;
            }
        }

        double plusleft1 = pt[0] = sum_x / sum_gx;
        double plusleft2 = pt[1] = sum_y / sum_gy;
        double plusright1 = pt[2] = ta0 + ta1*pt[0] + ta2*pt[1];
        double plusright2 = pt[3] = tb0 + tb1*pt[0] + tb2*pt[1];


        //if (re[1] > 0.5)
        //cout << "pt:" << pt[0] << "\t" << pt[1] << "\t" << pt[2] << "\t" << pt[3] << endl;
        pt[0] += P0[0];
        pt[1] += P0[1];
        pt[2] += P0[2];
        pt[3] += P0[3];
        pt[4] = re[0];



        //if (re[1] > 0.5)
        //cout << re[1] << "\t" << pt[0] << "\t" << pt[1] << endl;
        //cout << "matches point:" << pt[0] << "\t" << pt[1] << "\t" << pt[2] << "\t" << pt[3] << endl;
        if (re[1] > 0.95)
        {
            ncount_count_valid++;
            _left_lsm.push_back(Point2d(pt[0], pt[1]));
            _right_lsm.push_back(Point2d(pt[2], pt[3]));
            _hfactor.push_back(Point2d(th0, th1));
            _leftmove.push_back(Point2d(plusleft1, plusleft2));
            _rightmove.push_back(Point2d(plusright1, plusright2));
            _correvaluearray.push_back(re[1]);

            MatchResult _matchresult;
            _matchresult.corre = re[1];
            _matchresult._lx = pt[0];  _matchresult._ly = pt[1];
            _matchresult._rx = pt[2];  _matchresult._ry = pt[3];
            _matchresult._lpx = plusleft1; _matchresult._lpy = plusleft2;
            _matchresult._rpx = plusright1; _matchresult._rpy = plusright2;
            _lsmresult.push_back(_matchresult);

            double rateleft = 0.9 / sqrt(pow(plusleft1, 2.0) + pow(plusleft2, 2.0));
            _leftnewkey.push_back(Point2d(pt[0] + rateleft * plusleft1, pt[1] + rateleft * plusleft2));
            double rateright = 0.9 / sqrt(pow(plusright1, 2.0) + pow(plusright2, 2.0));
            _rightnewkey.push_back(Point2d(pt[2] + rateright * plusright1, pt[3] + rateright * plusright2));
        }

        delete my_kernel;
        delete[] coeffarray;// = new double[realwidth_ * realwidth_ * 8];
        delete[] staticarray;// = new double[realwidth_ * realheight_]; // 观测值
        delete[] tarray;// = new double[realwidth_ * realheight_ * 8];  // 常量矩阵
        delete[] pimgleft;
    }
}

bool ImgFacotryAPI::Mosaic_RigurouModel_Edge(MatchToolType& _lefttool, MatchToolType& _right, VirtualMachine* hei, int index_id)
{
	CSatOrbit helper;
	_lefttool._typeofFac = Enum_Match_RigurousModel;
	_right._typeofFac = Enum_Match_RigurousModel;
	SelfBuild(_lefttool, hei);
	SelfBuild(_right, hei);

	Point3d range_left[4] = { Point3d(379, 0, 0), Point3d(479, 0, 0), Point3d(379, 10785, 0), Point3d(479, 10785, 0) };
	Point3d range_right[4] = { Point3d(0, 0, 0), Point3d(100, 0, 0), Point3d(0, 10785, 0), Point3d(100, 10785, 0) };

	///获得左右影像兴趣区域对应的地面范围
	_lefttool.Translate(true, range_left, 4, index_id, false);
	_lefttool.Translate(true, range_right, 4, index_id+1, false);
	for (int i = 0; i < 4; i++)
	{
		cout << "left point:" << endl;
		cout << range_left[i] << endl;
		cout << "right point:" << endl;
		cout << range_right[i] << endl;
	}

	double lon_info[8];
	double lat_info[8];
	for (int i = 0; i < 4; i++)
	{
		lon_info[i * 2 + 0] = range_left[i].x;
		lon_info[i * 2 + 1] = range_right[i].x;

		lat_info[i * 2 + 0] = range_left[i].y;
		lat_info[i * 2 + 1] = range_right[i].y;
	}

	///找到影像经纬度范围
	double lon_min, lon_max;
	double lat_min, lat_max;
	safestring::Array_MinMax(lon_info, 8, lon_max, lon_min);
	safestring::Array_MinMax(lat_info, 8, lat_max, lat_min);

	cout << "经度范围:" << lon_min << "->" << lon_max << endl;
	cout << "纬度范围:" << lat_min << "->" << lat_max << endl;

	///确定分辨率
	double gsd_value[4] = { 0 };
	gsd_value[0] = static_cast<VirtualMachine*>(_lefttool._ptool)->distanceOnePixel(index_id, 379, 10, 100, 100, 100);
	gsd_value[1] = static_cast<VirtualMachine*>(_lefttool._ptool)->distanceOnePixel(index_id, 379, 10675, 100, 100, 100);
	gsd_value[2] = static_cast<VirtualMachine*>(_lefttool._ptool)->distanceOnePixel(index_id + 1, 0, 10, 100, 100, 100);
	gsd_value[3] = static_cast<VirtualMachine*>(_lefttool._ptool)->distanceOnePixel(index_id + 1, 0, 10675, 100, 100, 100);

	double gsdvalue = (gsd_value[0] + gsd_value[1] + gsd_value[2] + gsd_value[3]) / 4;


	///确定影像大小
	double x_range_value = (lat_max - lat_min) / gsdvalue;
	double y_range_value = (lon_max - lon_min) / gsdvalue;

	int x_range = x_range_value;
	int y_range = y_range_value;
	cout << "影像宽度:" << x_range << endl;
	cout << "影像宽度:" << y_range << endl;
	int pause_mode = 0;
	cout << gsd_value[0] << "\t" << gsd_value[1] << "\t" << gsd_value[2] << "\t" << gsd_value[3] << endl;
	cin >> pause_mode;

	GDALDriver* poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
	_right.outputpath = _lefttool.filepath + "_Mosaic.tiff";
	GDALDataset* poDataset = poDriver->Create(_right.outputpath.c_str(), x_range, y_range, 1, GDT_Float64, NULL);
	double* dataoutput = new double[x_range * y_range];
	memset(dataoutput, 0, sizeof(double)* x_range * y_range);
	GDALRasterBand* poBand = poDataset->GetRasterBand(1);

	Mat mat_left_ori = imread(_lefttool.filepath, CV_LOAD_IMAGE_ANYDEPTH);
	cout << _lefttool.filepath.c_str() << endl;
	Mat mat_left; mat_left_ori.convertTo(mat_left, CV_64F);
	cout << _right.filepath.c_str() << endl;
	int _hei_pause;
	cin >> _hei_pause;
	Mat mat_right_ori = imread(_right.filepath, CV_LOAD_IMAGE_ANYDEPTH);
	Mat mat_right; mat_right_ori.convertTo(mat_right, CV_64F);

	int ncount_range = x_range  * y_range;
	int ncounter_range = 0;
	Point3d* point_info = new Point3d[ncount_range];
	for (int j = 0; j < y_range; j++)
	{
		for (int i = 0; i < x_range; i++, ncounter_range++)
		{
			point_info[ncounter_range].y = lat_min + i * gsdvalue;
			point_info[ncounter_range].x = lon_min + j * gsdvalue;
			point_info[ncounter_range].z = 10;
		}
	}

	_lefttool.Translate(false, point_info, ncount_range, index_id, true);

	ResampleKernels* my_kernel = new ResampleKernels;
	my_kernel->Init(Raised_Cosine_6P);
	for (int i = 0; i < ncount_range; i++)
	{
		if (point_info[i].z == -9999)
			continue;
		else if (point_info[i].z == 0)
		{

			dataoutput[i] =  safeopencv::Resample(mat_left, Point2f(static_cast<float>(point_info[i].x), static_cast<float>(point_info[i].y)), my_kernel);

		}
		else if (point_info[i].z == 1)
		{
			dataoutput[i] = safeopencv::Resample(mat_right, Point2f(static_cast<float>(point_info[i].x), static_cast<float>(point_info[i].y)), my_kernel);
		}
	}

	_lefttool.Translate(false, point_info, ncount_range, index_id, true);

	poBand->RasterIO(GF_Write, 0, 0, x_range, y_range, dataoutput, x_range, y_range, GDT_Float64, 0, 0);
	GDALClose(poDataset);
	delete[] point_info;
	delete my_kernel;
	return true;
}


bool ImgFacotryAPI::LSM_Match(MatchToolType leftinfi, MatchToolType rightinfi, int windowsize)
{

    vector<MatchResult> _resultcontainer;
    vector<Point2d> leftsecond; 
    vector<Point2d> rightsecond;

    LSM_Match(leftinfi, rightinfi, _leftkey_, _rightkey_, leftsecond, rightsecond, _resultcontainer, windowsize);
    if (true)
    {
        vector<Point2d> leftthird;
        vector<Point2d> rightthird;


        LSM_Match(leftinfi, rightinfi, leftsecond, rightsecond, leftthird, rightthird, _resultcontainer, windowsize);
        vector<Point2d> leftnext;
        vector<Point2d> rightnext;

        LSM_Match(leftinfi, rightinfi, leftthird, rightthird, leftnext, rightnext, _resultcontainer, windowsize);

    }
    cout << "match point number" << endl;

    ofstream ofsmatch("c:\\Temp\\match\\heimatch.xls");
    cout << _resultcontainer.size() << endl;
    ofsmatch.precision(20);
    vector<Point2d> yarrayclib;
    for (int i = 0; i < _resultcontainer.size(); i++)
    {
        ofsmatch << _resultcontainer[i]._lx << "\t" << _resultcontainer[i]._ly << "\t" << _resultcontainer[i]._rx << "\t" << _resultcontainer[i]._ry << "\t" << "\t" << _resultcontainer[i]._lpx << "\t" << _resultcontainer[i]._lpy << "\t " << _resultcontainer[i]._rpx << "\t" << _resultcontainer[i]._rpy << "\t" << _resultcontainer[i].corre
            << "\t" << _resultcontainer[i]._rx - _resultcontainer[i]._lx << "\t" << _resultcontainer[i]._ry - _resultcontainer[i]._ly << endl;
        _heitemp.push_back(Point2d(_resultcontainer[i]._ly, _resultcontainer[i]._ly - _resultcontainer[i]._ry));
        _heitemp_x.push_back(Point2d(_resultcontainer[i]._lx, _resultcontainer[i]._ly - _resultcontainer[i]._ry));
    }

    return true;
}

bool ImgFacotryAPI::Template_PreMatch(string match_path)
{
    //清空原有存储
    _leftkey_.clear();
    _rightkey_.clear();
    ifstream ifs(match_path.c_str());

    if (!ifs.is_open())
    {
        cout << "精匹配结果文件不存在，请重新计算或者使用特征匹配粗匹配结果" << endl;
        return false;
    }

    else
    {
        cout << "读取匹配控制点...." << endl;
        do
        {
//210.9727807195454	617.65362414936408	208.9727807195454	617.65362414936408		-0.02721928045460284	-0.34637585063591253	 -0.02721928045460284	-0.34637585063591253	0.97769335571224025
            char reader[1024];
            ifs.getline(reader, 1024);
            double a[11];
            sscanf_s(reader, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &a[0], &a[1], &a[2], &a[3], &a[4], &a[5], &a[6], &a[7], &a[8],&a[9],&a[10]);
            _leftkey_.push_back(Point2d(a[0], a[1]));
            _rightkey_.push_back(Point2d(a[2], a[3]));
        } while (ifs.peek()!=EOF);

        cout << "匹配控制点读取完毕" << endl;

        ///point interpolation
        cout << "开始拟合匹配曲线..." << endl;


        safeopencv::Match_Curve(_leftkey_, _rightkey_, "nonesense", "nonesense",x2xinfo,y2yinfo);

        cout << "匹配曲线拟合完毕" << endl;
        
        ///根据匹配结果输入模版匹配初值



    }
    return false;
}


bool ImgFacotryAPI::PreMatch(MatchToolType& _lefttool, MatchToolType& _right, VirtualMachine* _leftModel, VirtualMachine* _rightModel)
{
    CSatOrbit helper;
    SelfBuild(_lefttool, _leftModel);
    SelfBuild(_right, _rightModel);





    cv::Rect hei;
    Point2f srcpt[3];
    Point2f dstpt[3];
    if (_lefttool._typeofFac == Enum_Match_RigurousModel)
    {
        srcpt[0].x = 0; srcpt[0].y = 0;
        int leftwidth = static_cast<VirtualMachine*>(_lefttool._ptool)->SingWidth();
        int leftheight = static_cast<VirtualMachine*>(_lefttool._ptool)->SingHeight();
        double rectXInfo[4] = { 0, leftwidth - 1, 0, leftwidth-1 };
        double rectYInfo[4] = { 0, 0, leftheight -1, leftheight-1 };
        Translate(_lefttool, _right, 0, rectXInfo, rectYInfo, 4);
        for (int t = 0; t < 4; t++)
        {
            Point2d heishowit(rectXInfo[t], rectYInfo[t]);
            cout << heishowit << endl;
        }

        double x_Min, x_Max;
        double y_Min, y_Max;
        int sizearray = 4;
        safestring::Array_MinMax(rectXInfo, sizearray, x_Max, x_Min);
        safestring::Array_MinMax(rectYInfo, sizearray, y_Max, y_Min);

        GDALDataset* poRight = (GDALDataset*)GDALOpen(_right.filepath.c_str(), GA_ReadOnly);
        int rightwidth = poRight->GetRasterXSize();
        int rightheight = poRight->GetRasterYSize();
        ///确定需要搜寻的范围，减小内存消耗
        cout << y_Min << endl;
        x_Min = ((x_Min - 10) < 0) ? 0 : x_Min - 10;
        x_Max = ((x_Max + 10) > rightwidth - 1) ? rightwidth - 1 : x_Max + 10;
        y_Min = ((y_Min - 10) < 0) ? 0 : y_Min - 10;
        y_Max = ((y_Max + 10) > rightheight - 1) ? rightheight - 1 : y_Max + 10;


        //后面的就交给GDAL
        GDALDriver* poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
        _right.outputpath = _right.filepath + "PreMatch.tiff";
        GDALDataset* poDataset = poDriver->Create(_right.outputpath.c_str(), leftwidth, leftheight, 1, GDT_Float64, NULL);
        double* dataoutput = new double[leftwidth * leftheight];
        memset(dataoutput, 0, sizeof(double)* leftwidth * leftheight);
        GDALRasterBand* poBand = poDataset->GetRasterBand(1);
        GDALRasterBand* poBandOri = poRight->GetRasterBand(1);

        int block_Xsize, block_Ysize;
        block_Xsize = x_Max - x_Min + 1;
        block_Ysize = y_Max - y_Min + 1;
        int* ori_Block_info = new int[block_Xsize * block_Ysize];
        poBandOri->RasterIO(GF_Read, x_Min, y_Min, block_Xsize, block_Ysize, ori_Block_info, block_Xsize, block_Ysize, GDT_Int32, 0, 0);
        cout << "block_X_size:" << block_Xsize << endl;
        cout << "block_Y_szie:" << block_Ysize << endl;
        cout << "block_X_start:" << x_Min << endl;
        cout << "block_Y_start:" << y_Min << endl;
        //#pragma omp parallel for
        for (int j = 0; j < leftheight; j++)
        {
            cout << "Progress:" << j << "/" << leftheight << "\r";
            ///设定采样核
            if (j > 1000) break;
            ResampleKernels* my_kernel = new ResampleKernels;
            my_kernel->Init(Raised_Cosine_6P);
            int kernelsize = my_kernel->GetCurKernelSize();
            int* kerneldata = new int[kernelsize * kernelsize];
            int pixelsize = kernelsize * kernelsize;
            double* kernelinfo = new double[kernelsize * kernelsize];
            for (int i = 0; i < leftwidth; i++)
            {
                double x_info[1] = { i }, y_info[1] = { j };
                Translate(_lefttool, _right, 0, x_info, y_info, 1);

                ///将影像坐标转换为Block坐标
                double block_loc_x = x_info[0] - x_Min;
                double block_loc_y = y_info[0] - y_Min;
                int up_index_x = block_loc_x - kernelsize + 1;
                int up_index_y = block_loc_y - kernelsize + 1;
                int down_index_x = block_loc_x + kernelsize - 1;
                int down_index_y = block_loc_y + kernelsize - 1;



                ///尽量在for循环内部避免if语句，将边界判断放在外部
                ///采用方式，针对不同采样核大小，给影像镶嵌一圈对称的边界，
                for (int q = 0; q < kernelsize; q++)
                {
                    memcpy(kerneldata + q * kernelsize, ori_Block_info + (up_index_y + q)* block_Xsize + up_index_x, sizeof(int)*kernelsize);
                }

                my_kernel->GetCurPixelResampleKernels(block_loc_x, block_loc_y, kernelinfo);
                double sum_all_block_pixel = 0;
                for (int q = 0; q < pixelsize; q++)
                {
                    sum_all_block_pixel += kerneldata[q] * kernelinfo[q];
                }

                dataoutput[j * leftwidth + i] = sum_all_block_pixel * 1463.800049 / 26858;
                //cout << dataoutput[j * leftwidth + i] << endl;
            }

            delete[] kerneldata;
            delete[] kernelinfo;
        }

        poBand->RasterIO(GF_Write, 0, 0, leftwidth, leftheight, dataoutput, leftwidth, leftheight, GDT_Float64, 0, 0);
        delete[] dataoutput;
        GDALClose(poDataset);
        GDALClose(poRight);
    }


    return true;
}


ImgFacotryAPI::~ImgFacotryAPI()
{
    if (_imgcontent != nullptr)
    {
        delete[] _imgcontent;
        _imgcontent = nullptr;
    }
}

unsigned int* ImgFacotryAPI::IRow(unsigned int _rowindex)
{
    return _imgcontent + _rowindex * _width;
}

//bool ImgFacotryAPI::IOutPut(bool radiometrix)


void ImgFacotryAPI::Mosaic(vector<string> imagenames, int width, int height, bool align)
{
    GDALDataType dataType;
    if (_singlebyte <= 8)
    {
        dataType = GDT_Byte;
    }

    else if (_singlebyte > 8 && _singlebyte <= 16)
    {
        dataType = GDT_UInt16;
    }

    else if (_singlebyte > 16 && _singlebyte <= 32)
    {
        dataType = GDT_UInt32;
    }

    else if (_singlebyte > 32 && _singlebyte <= 64)
    {
        dataType = GDT_Float64;
    }

    else
    {
        return;
    }

    string pathof = "Whole.tif";
    GDALDriver* poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    GDALDataset* poDataset = poDriver->Create(pathof.c_str(), width, height, 1, dataType, NULL);
    GDALRasterBand* poBand = poDataset->GetRasterBand(1);
    int eachwidth, eachheight; eachwidth = eachheight = 0;
    if (!align)
    {
        eachheight = height;
        eachwidth = width / imagenames.size();
    }
    else
    {
        eachheight = height / imagenames.size();
        eachwidth = width;
    }
    for (int i = 0; i < imagenames.size(); i++)
    {
        GDALDataset* pgdaltemp = (GDALDataset*)GDALOpen(imagenames[i].c_str(), GA_ReadOnly);
        GDALRasterBand* pobandtemp = pgdaltemp->GetRasterBand(1);
        double* infocont = new double[eachwidth * eachheight];
        pobandtemp->RasterIO(GF_Read, 0, 0, eachwidth, eachheight, infocont, eachwidth, eachheight, dataType, 0, 0);

        if (!align)
        {
            poBand->RasterIO(GF_Write, eachwidth*i, 0, eachwidth, eachheight, infocont, eachwidth, eachheight, dataType, 0, 0);
        }
        else
        {
            poBand->RasterIO(GF_Write, 0, eachheight*i, eachwidth, eachheight, infocont, eachwidth, eachheight, dataType, 0, 0);

        }
        delete[] infocont;
        GDALClose(pgdaltemp);
    }
    GDALClose(poDataset);
}

bool ImgFacotryAPI::Phase_Correction(Mat matleft, Mat matright, double& maxvalue, Point2d& pointmove)
{
    if (matleft.cols != matright.cols)
    {
        cout << "错误:分块影像尺寸不一致!" << endl;
        return false;
    }

    if (matleft.rows != matright.rows)
    {
        cout << "错误:影像尺寸不一致!" << endl;
        return false;
    }

    if (matleft.cols != matleft.rows)
    {
        cout << "请保持窗口宽高一致" << endl;
        return false;
    }

    int imgwidth, imgheight;
    imgwidth = imgheight = matleft.cols;

    void* pvoidleft = static_cast<void*>(matleft.data);
    double* pimgleft = static_cast<double*>(pvoidleft);
    void* pvoidright = static_cast<void*>(matright.data);
    double* pimgright = static_cast<double*>(pvoidright);

    fftw_complex* complex_left = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)* imgwidth * imgheight));
    fftw_complex* complex_right = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)* imgwidth * imgheight));
    fftw_complex* complex_res = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)* imgwidth * imgheight));

    fftw_plan fft_plan_left = fftw_plan_dft_1d(imgwidth * imgheight, complex_left, complex_left, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan fft_plan_right = fftw_plan_dft_1d(imgwidth * imgheight, complex_right, complex_right, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan fft_plan_res = fftw_plan_dft_1d(imgwidth * imgheight, complex_res, complex_res, FFTW_BACKWARD, FFTW_ESTIMATE);


    int k = 0;
    for (int j = 0; j < imgheight; j++)
    {
        for (int i = 0; i < imgwidth; i++, k++)
        {
            complex_left[k][0] = pimgleft[k];
            complex_left[k][1] = 0;

            complex_right[k][0] = pimgright[k];
            complex_right[k][1] = 0;
        }
    }

    fftw_execute(fft_plan_left);
    fftw_execute(fft_plan_right);

    int fft_size = imgwidth *imgheight;
    for (int i = 0; i < fft_size; i++) {
        complex_res[i][0] = (complex_right[i][0] * complex_left[i][0]) - (complex_right[i][1] * (-complex_left[i][1]));
        complex_res[i][1] = (complex_right[i][0] * (-complex_left[i][1])) + (complex_right[i][1] * complex_left[i][0]);

        double tmp = sqrt(pow(complex_res[i][0], 2.0) + pow(complex_res[i][1], 2.0));

        complex_res[i][0] /= tmp;
        complex_res[i][1] /= tmp;

    }

    fftw_execute(fft_plan_res);
    Mat fftwresult(imgwidth, imgheight, CV_32F);

    void* pvoidresult = static_cast<void*>(fftwresult.data);
    float* pImgresult = static_cast<float*>(pvoidresult);
    /* normalize and copy to result image */
    double calc_max = -100;
    int index_i = 0;
    int index_j = 0;
    int index_min_i = 0;
    int index_min_j = 0;
    double calc_min = 100;
    for (int i = 0; i < fft_size; i++) {
        pImgresult[i] = complex_res[i][0] / (double)fft_size;
        if (pImgresult[i] > calc_max)
        {
            index_j = i / imgwidth;
            index_i = i % imgwidth;
        }
        if (pImgresult[i] < calc_min)
        {
            index_min_i = i / imgwidth;
            index_min_j = i % imgwidth;
            calc_min = pImgresult[i];
        }
        calc_max = (pImgresult[i] > calc_max) ? pImgresult[i] : calc_max;

    }


    maxvalue = calc_max;
    pointmove.x = (index_i);
    pointmove.y = (index_j);

    Point2d min_max_Loc[2];
    Point min_max_ori[2];
    min_max_ori[0].x = index_i; min_max_ori[0].y = index_j;
    min_max_ori[1].x = index_min_i; min_max_ori[1].y = index_min_j;
    safeopencv::MinMaxloc_subpixel(min_max_Loc, fftwresult, min_max_ori, 0);
    pointmove.x = min_max_Loc[0].x; pointmove.y = min_max_Loc[0].y;


    if (min_max_Loc[0].x > imgwidth / 2 || min_max_Loc[0].y > imgheight / 2)
    {
        if (!is_shifted)
        {
            is_shifted = !is_shifted;
            Phase_Correction(matright, matleft, maxvalue, pointmove);
            pointmove.x = -pointmove.x;
            pointmove.y = -pointmove.y;
        }
    }

    /*
    char readernow[1024];
    heivaluecount++;
    sprintf_s(readernow,"c:\\Temp\\fftw_%d.tiff", heivaluecount);
    GDALDriver* poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    GDALDataset* poDatasetNew = poDriver->Create(readernow, imgwidth, imgheight, 1, GDT_Float32, NULL);
    GDALRasterBand* poBandNew = poDatasetNew->GetRasterBand(1);
    poBandNew->RasterIO(GF_Write, 0, 0, imgwidth, imgheight, pImgresult, imgwidth, imgheight, GDT_Float32, 0, 0);



    delete poDatasetNew;
    */
    return true;
}

bool ImgFacotryAPI::Phase_Correction(string _leftpath, string _rightpath)
{
    GDALDataset* poDataleft = (GDALDataset*)GDALOpen(_leftpath.c_str(), GA_ReadOnly);
    GDALDataset* poDataright = (GDALDataset*)GDALOpen(_rightpath.c_str(), GA_ReadOnly);
    GDALRasterBand* poBandleft = poDataleft->GetRasterBand(1);
    GDALRasterBand* poBandright = poDataright->GetRasterBand(1);

    if (poDataleft->GetRasterXSize() != poDataright->GetRasterXSize())
    {
        cout << "错误:影像尺寸不一致!" << endl;
        return false;
    }

    if (poDataleft->GetRasterYSize() != poDataright->GetRasterYSize())
    {
        cout << "错误:影像尺寸不一致!" << endl;
        return false;
    }
    int imgwidth = poDataleft->GetRasterXSize();
    int imgheight = poDataleft->GetRasterYSize();

    double* leftdata = new double[imgwidth * imgheight];
    int* rightdataori = new int[imgwidth * imgheight];
    double* rightdata = new double[imgwidth * imgheight];

    poBandleft->RasterIO(GF_Read, 0, 0, imgwidth, imgheight, leftdata, imgwidth, imgheight, GDT_Float64, 1, 1);
    poBandright->RasterIO(GF_Read, 0, 0, imgwidth, imgheight, rightdataori, imgwidth, imgheight, GDT_Int32, 1, 1);
    for (int j = 0; j < imgheight; j++)
    for (int i = 0; i < imgwidth; i++)
    {
        rightdata[j * imgwidth + i] = rightdataori[j*imgwidth + i];
    }
    delete[] rightdataori;

    fftw_complex* complex_left = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)* imgwidth * imgheight));
    fftw_complex* complex_right = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)* imgwidth * imgheight));
    fftw_complex* complex_res = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)* imgwidth * imgheight));

    fftw_plan fft_plan_left = fftw_plan_dft_1d(imgwidth * imgheight, complex_left, complex_left, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan fft_plan_right = fftw_plan_dft_1d(imgwidth * imgheight, complex_right, complex_right, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan fft_plan_res = fftw_plan_dft_1d(imgwidth * imgheight, complex_res, complex_res, FFTW_BACKWARD, FFTW_ESTIMATE);

    ///input data
    int k = 0;
    for (int j = 0; j < imgheight; j++)
    {
        for (int i = 0; i < imgwidth; i++, k++)
        {
            complex_left[k][0] = leftdata[k];
            complex_left[k][1] = 0;

            complex_right[k][0] = rightdata[k];
            complex_right[k][1] = 0;
        }
    }

    fftw_execute(fft_plan_left);
    fftw_execute(fft_plan_right);

    int fft_size = imgwidth *imgheight;
    for (int i = 0; i < fft_size; i++) {
        complex_res[i][0] = (complex_right[i][0] * complex_left[i][0]) - (complex_right[i][1] * (-complex_left[i][1]));
        complex_res[i][1] = (complex_right[i][0] * (-complex_left[i][1])) + (complex_right[i][1] * complex_left[i][0]);

        double tmp = sqrt(pow(complex_res[i][0], 2.0) + pow(complex_res[i][1], 2.0));

        complex_res[i][0] /= tmp;
        complex_res[i][1] /= tmp;
    }

    fftw_execute(fft_plan_res);

    string outputpath = _leftpath + "phase_correlation.tiff";
    GDALDriver* poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    GDALDataset* poDatasetNew = poDriver->Create(outputpath.c_str(), imgwidth, imgheight, 1, GDT_Float64, NULL);
    GDALRasterBand* poBandNew = poDatasetNew->GetRasterBand(1);
    double * datainfo = new double[imgwidth * imgheight];
    /* normalize and copy to result image */
    for (int i = 0; i < fft_size; i++) {
        datainfo[i] = complex_res[i][0] / (double)fft_size;
    }

    poBandNew->RasterIO(GF_Write, 0, 0, imgwidth, imgheight, datainfo, imgwidth, imgheight, GDT_Float64, 1, 1);

    GDALClose(poDataleft);
    GDALClose(poDataright);

    fftw_destroy_plan(fft_plan_left);
    fftw_destroy_plan(fft_plan_right);
    fftw_destroy_plan(fft_plan_res);
    fftw_free(complex_left);
    fftw_free(complex_right);
    fftw_free(complex_res);

    delete[] leftdata;
    delete[] rightdata;
    delete[] datainfo;
}

bool Match_Phase_Correction(MatchToolType& _lefttool, MatchToolType& _righttool, int windowsize, int windowmove, double thresold)
{
    int xstep = static_cast<VirtualMachine*>(_lefttool._ptool)->SingWidth();
    int ystep = static_cast<VirtualMachine*>(_lefttool._ptool)->SingHeight();
    int halfsize = windowsize / 2;
    int xedge = xstep - windowsize;
    int yedge = ystep - windowsize;
    int xmovestep = xedge / windowmove;
    int ymovestep = yedge / windowmove;

    ///保持影像原有的级数不变
    Mat imgleft = imread(_lefttool.outputpath, CV_LOAD_IMAGE_ANYDEPTH);
    Mat imgright = imread(_righttool.outputpath, CV_LOAD_IMAGE_ANYDEPTH);

    double transfromleft2right[6];
    for (int j = 0; j < ymovestep; j++)
    {
        for (int i = 0; i < xmovestep; i++)
        {
            //获取分块 
            Rect_<int> leftrange(i * windowsize, j * windowsize, windowsize, windowsize);
            int leftImg_Top_X = i * windowsize; int rightImg_Top_y = j * windowsize;
            int rightImg_Top_X, rightImg_Top_Y;
            rightImg_Top_X = leftImg_Top_X * transfromleft2right[0] + rightImg_Top_y * transfromleft2right[1] + transfromleft2right[2];
            rightImg_Top_Y = leftImg_Top_X * transfromleft2right[3] + rightImg_Top_y * transfromleft2right[4] + transfromleft2right[5];
            Rect_<int> rightrange(rightImg_Top_X, rightImg_Top_Y, windowsize, windowsize);
            Mat arealeft = imgleft(leftrange);
            Mat arearight = imgright(rightrange);



        }
    }

    return true;
}

bool  ImgFacotryAPI::StitchFile(vector<string> iamgenames, string output, int eachwidht, int eachheight)
{
    vector<Mat> affineinfo;
    float affinepara[3][2] = { 0 };
    Mat mataffine;
    cout << mataffine << endl;
    MatchAndCompare("C:\\Temp\\", &mataffine);
    for (int i = 0; i < iamgenames.size() - 1; i++)
    {
        /*
        cv::Ptr<cv::Mat> affine;// = cv::Mat::create(2, 3, CV_32F);
        Rect rectleft(eachwidht - 120, 1000, 99, 5000);
        Rect rectright(0, 1000, 99, 5000);
        MatchAndCompare(iamgenames[i], iamgenames[i + 1], affine,rectleft,rectright);*/
        affineinfo.push_back(mataffine);
    }

    vector<Point2f> stichpoints;
    Point2f hei1(0, 0);
    stichpoints.push_back(hei1);
    for (int i = 1; i < iamgenames.size(); i++)
    {
        Point2f hei(1, 1);
        for (int j = i - 1; j >= 0; j--)
        {
            float heitemp[2] = { hei.x, hei.y };
            //Mat matemp(2, 1, CV_32F, &heitemp);
            float heiresult[2] = { 0 };
            //cout << affineinfo[j] << endl;
            double a[2][3] = { 0 };
            for (int t = 0; t < 2; t++)
            {
                for (int q = 0; q < 3; q++)
                {
                    a[t][q] = affineinfo[j].at<double>(t, q);
                }
            }
            hei.x = hei.x * a[0][0] + hei.y * a[0][1] + a[0][2];
            hei.y = hei.x * a[1][0] + hei.y * a[1][1] + a[1][2];
        }

        stichpoints.push_back(hei);
    }

    Mosaic<unsigned short>(iamgenames, stichpoints, eachwidht, eachheight);
    return true;
}

void ImgFacotryAPI::SplitByCol(string imagename, int singleccdbyte)
{
    GDALDataset* poDataSet = (GDALDataset*)GDALOpen(imagename.c_str(), GA_ReadOnly);
    GDALRasterBand* poBand = poDataSet->GetRasterBand(1);
    int sizeofcol = poDataSet->GetRasterXSize() / singleccdbyte;
    GDALDriver* poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    int heightfo = poDataSet->GetRasterYSize();
    int widthfo = poDataSet->GetRasterXSize();
    vector<unsigned short*> oriinfo;
    unsigned short* _datainfoori = new unsigned short[heightfo * widthfo];
    poBand->RasterIO(GF_Read, 0, 0, widthfo, heightfo, _datainfoori, widthfo, heightfo, GDT_UInt16, 0, 0);
    oriinfo.resize(widthfo);


    for (int i = 0; i < widthfo; i++)
    {
        unsigned short* _datainfo = new unsigned short[heightfo];
        //poBand->RasterIO(GF_Read, i, 0, 1, heightfo, _datainfo, 1, heightfo, GDT_UInt16, 0, 0);
        for (int j = 0; j < heightfo; j++)
        {
            _datainfo[j] = _datainfoori[j * widthfo + i];
        }
        oriinfo[i] = _datainfo;
    }
    delete[] _datainfoori;


    for (int i = 0; i < singleccdbyte; i++)
    {
        char reader[1024];
        sprintf_s(reader, "%d_new.tif", i);
        GDALDataset* poDatasetNew = poDriver->Create((imagename + static_cast<string>(reader)).c_str(), heightfo, sizeofcol, 1, GDT_UInt16, NULL);
        GDALRasterBand* poBandNew = poDatasetNew->GetRasterBand(1);
        unsigned short * datainfo = new unsigned short[sizeofcol * heightfo];
        for (int j = 0; j < sizeofcol; j++)
        {
            memcpy(datainfo + heightfo * j, oriinfo[j * singleccdbyte + i], sizeof(unsigned short)*heightfo);
        }

        poBandNew->RasterIO(GF_Write, 0, 0, heightfo, sizeofcol, datainfo, heightfo, sizeofcol, GDT_UInt16, 0, 0);
        delete[]datainfo;
        GDALClose(poDatasetNew);
    }

    GDALClose(poDataSet);

    for (int i = 0; i < oriinfo.size(); i++)
    {
        if (oriinfo[i] != nullptr)
        {
            delete[] oriinfo[i];
            oriinfo[i] = nullptr;
        }
    }
}


#include "safestring.h"
//匹配影像，给出仿射系数矩阵
void ImgFacotryAPI::MatchAndCompare(string rootdir, Mat* heimat)
{
    //327.0000 46.0000 204.1558 27.4547
    vector<string> mpointre;
    vector<Point2f> leftpoints;
    vector<Point2f> rightpoint;
    safestring::FindFileSInDir(rootdir, mpointre, "*.mpoints");
    for (int i = 0; i < mpointre.size(); i++)
    {
        ifstream ifs;
        ifs.open(mpointre[i].c_str());
        char reader[1024];
        ifs.getline(reader, 1024);
        int pointsize = 0;
        sscanf_s(reader, "%d", &pointsize);
        for (int j = 0; j < pointsize; j++)
        {
            ifs.getline(reader, 1024);
            float x1, y1, x2, y2;
            sscanf_s(reader, "%f %f %f %f", &x1, &y1, &x2, &y2);
            if (x2 > 22.0)
                continue;
            if (y2 < y1)
                continue;
            Point2f p1, p2;
            p1.x = x1; p1.y = y1;
            p2.x = x2; p2.y = y2;
            leftpoints.push_back(p1);
            rightpoint.push_back(p2);
        }
    }
    double heivaluedouble[3][2] = { 0 };
    //Mat hei(3, 2, CV_64F, &heivaluedouble);
    Mat hei = cv::estimateRigidTransform(rightpoint, leftpoints, false);
    cout << hei << endl;
    hei.copyTo(*heimat);
}

void ImgFacotryAPI::MatchAndCompare(string _left, string _right, cv::Ptr<cv::Mat> matinfo, cv::Rect rectleft, cv::Rect rectRight)
{
    initModule_nonfree();
    cv::Ptr<cv::FeatureDetector> detector = cv::FeatureDetector::create("HARRIS");
    //printParams(detector);

    //detector->set("edgeThreshold", 1);
    cv::Ptr<cv::DescriptorExtractor> descriptor = cv::DescriptorExtractor::create("SIFT");

    //detector->
    ///获得关键点
    std::vector<cv::KeyPoint> keypoints1, keypoints2, keypionts3;

    Mat img_1 = imread(_left.c_str(), CV_LOAD_IMAGE_ANYDEPTH);
    Mat img_2 = imread(_right.c_str(), CV_LOAD_IMAGE_ANYDEPTH);

    //imshow("hei",picu8bit);
    //waitKey(0);
    cv::initModule_nonfree();
    //detector->detect(picu8bit(rectleft), keypionts3);

    if (!img_1.data || !img_1.data)
    {
        cout << "影像读取出错" << endl;
        return;
    }

    cv::initModule_nonfree();
    //detector->detect(img_1, keypoints1);
    detector->detect(img_1(rectleft), keypoints1);
    detector->detect(img_2(rectRight), keypoints2);

    ///extract features
    Mat descriptions_1, descriptions_2;

    descriptor->compute(img_1, keypoints1, descriptions_1);
    descriptor->compute(img_2, keypoints2, descriptions_2);

    BFMatcher matcher;
    vector<DMatch> matches1;
    matcher.match(descriptions_1, descriptions_2, matches1);

    Mat img_matches;
    drawMatches(img_1, keypoints1, img_2, keypoints2,
        matches1, img_matches, Scalar::all(-1), Scalar::all(-1),
        vector<char>(), DrawMatchesFlags::DRAW_RICH_KEYPOINTS);
    imwrite("d:\\matchresult.tif", img_matches);
    vector<Point2f> _leftpoints;
    vector<Point2f> _rightpoints;

    for (int i = 0; i < matches1.size(); i++)
    {
        _leftpoints.push_back(keypoints1[matches1[i].queryIdx].pt);
        _rightpoints.push_back(keypoints2[matches1[i].trainIdx].pt);
    }

    Mat hei = cv::estimateRigidTransform(_leftpoints, _rightpoints, false);
    hei.copyTo(*matinfo);
}



///进行直方图匹配，获得查找表
void ImgFacotryAPI::HistoMatching(Mat* info, MatND hisinfo, int singlebytenumber, unsigned short* T)
{
    int bins = static_cast<int>(pow(2, singlebytenumber));
    int hist_size[] = { bins };
    float range[] = { 0, static_cast<float>(bins) };
    const float* ranges[] = { range };
    Mat heiinfo;

    MatND hiswhole; int channels[] = { 0 };
    /* calcHist( &gray, 1, channels, Mat(), // do not use mask
    hist, 1, hist_size, ranges,
    true, // the histogram is uniform
    false );  */
    //int hei = info->at<unsigned short>(1000,5);
    cvtColor(*info, heiinfo, CV_BGR2GRAY);
    imwrite("hei.tiff", heiinfo);
    calcHist(&heiinfo, 1, channels, Mat(), hiswhole, 1, hist_size, ranges, true, false);

    long val_1 = 0.0;
    long val_2 = 0.0;
    //unsigned short* T = new unsigned short[bins];
    long* S = new long[bins];
    long* G = new long[bins];
    for (int index = 0; index < bins; index++)
    {
        val_1 += hiswhole.at<float>(index);
        val_2 += hisinfo.at<float>(index);
        S[index] = val_1;
        G[index] = val_2;
    }

    ///归一化
    for (int index = 0; index < bins; index++)
    {
        S[index] = S[index] / S[bins - 1];
        G[index] = G[index] / G[bins - 1];
    }

    double min_val = bins;
    int PG = 0;
    for (int i = 0; i < bins; i++)
    {
        min_val = 1;
        for (int j = 0; j < bins; j++)
        {
            if (fabs((double)(G[j] - S[i])) < min_val && (G[j] - S[i]) >= 0)
            {
                min_val = fabs(double((G[j] - S[i])));
                PG = j;
            }
        }

        T[i] = PG;
    }

    delete[] S;
    delete[] G;
}

void ImgFacotryAPI::DrawHist(const MatND& b_hist)
{
    // Draw the histograms for B, G and R
    int hist_w = 4096; int hist_h = 400;
    int bin_w = cvRound((double)hist_w / 4096);

    Mat histImage(hist_h, hist_w, CV_8U);

    /// Normalize the result to [ 0, histImage.rows ]
    normalize(b_hist, b_hist, 0, histImage.rows, NORM_MINMAX, -1, Mat());


    /// Draw for each channel
    for (int i = 1; i < 4096; i++)
    {
        line(histImage, Point(bin_w*(i - 1), hist_h - cvRound(b_hist.at<float>(i - 1))),
            Point(bin_w*(i), hist_h - cvRound(b_hist.at<float>(i))),
            Scalar(255, 0, 0), 2, 8, 0);

    }

    imwrite("d:\\match_hist.bmp", histImage);

}

void ImgFacotryAPI::HistoMatching(string matinfo, MatND hisinfo, int singlebytenumber, float* T)
{
    int bins = static_cast<int>(pow(2, singlebytenumber));
    int hist_size[] = { bins };
    float range[] = { 0, static_cast<float>(bins) };
    const float* ranges[] = { range };
    Mat heiinfo_ori = imread(matinfo, CV_LOAD_IMAGE_ANYDEPTH);
    Mat heiinfo; heiinfo_ori.convertTo(heiinfo, CV_16U);
    MatND hiswhole; int channels[] = { 0 };
    /* calcHist( &gray, 1, channels, Mat(), // do not use mask
    hist, 1, hist_size, ranges,
    true, // the histogram is uniform
    false );  */
    //int hei = info->at<unsigned short>(1000,5);

    calcHist(&heiinfo, 1, channels, Mat(), hiswhole, 1, hist_size, ranges, true, false);




    long double val_1 = 0.0;
    long double val_2 = 0.0;
    //unsigned short* T = new unsigned short[bins];
    long double* S = new long double[bins];
    long double* G = new long double[bins];
    for (int index = 1; index < bins; index++)
    {
        val_1 += hiswhole.at<float>(index);
        val_2 += hisinfo.at<float>(index);
        S[index] = val_1;
        G[index] = val_2;
    }

    S[0] = 0;
    G[0] = 0;
    for (int index = 0; index < bins; index++)
    {
        S[index] = S[index] / S[bins - 1];
        G[index] = G[index] / G[bins - 1];
    }


    int PG = 0;
    int sum = 0;
    ofstream fs1("c:\\Temp\\match\\3.txt");
    for (int i = 0; i < bins; i++)
    {
        int indexnow = 0;
        vector<Point2d> cureinfo;
        cureinfo.resize(bins);
        ///本来就是0就没拯救的必要了
        if (S[i] == 0 || i == 0)
        {
            T[i] = 0;
            continue;
        }
        for (int j = 0; j < bins; j++)
        {
            cureinfo[j].x = j;
            cureinfo[j].y = fabs(S[i] - G[j]);
        }

        double min_val, max_val, min_loc, max_loc;
        safeopencv::Curve_MinMax(cureinfo, max_val, min_val, max_loc, min_loc);
        min_loc = (min_loc > bins - 1) ? bins - 1 : min_loc;
        min_loc = (min_loc < 0) ? 0 : min_loc;

        T[i] = min_loc;

        fs1 << T[i] << endl;
    }

    fs1.close();
    delete[] S;
    delete[] G;
}

void ImgFacotryAPI::ImageScale(string path, double topx, double topy, double width, double height, double scalex, double  scaley)
{
    GDALDataset* poDataset = (GDALDataset*)GDALOpen(path.c_str(), GA_ReadOnly);
    GDALRasterBand* poBand = poDataset->GetRasterBand(1);
    int imagewidth = width * scalex;
    int imageheight = height * scaley;

    unsigned short* datainfo = new unsigned short[imagewidth * imageheight];
    memset(datainfo, 0, sizeof(unsigned short)*  imagewidth * imageheight);
    poBand->RasterIO(GF_Read, topx, topy, width, height, datainfo, imagewidth, imageheight, GDT_UInt16, 0, 0);

    GDALDriver* poDriver = GetGDALDriverManager()->GetDriverByName("ENVI");
    GDALDataset* poDatasetNew = poDriver->Create((path + "hei.img").c_str(), imagewidth, imageheight, 1, GDT_UInt16, NULL);
    GDALRasterBand* poBandNew = poDatasetNew->GetRasterBand(1);

    poBandNew->RasterIO(GF_Write, 0, 0, imagewidth, imageheight, datainfo, imagewidth, imageheight, GDT_UInt16, 0, 0);

    GDALClose(poDataset);
    GDALClose(poDatasetNew);
}

///Deprecated code
/*
CSatOrbit helper;
if (windowsize % 2 == 0)
{
cout << "请将最小二乘窗口大小设置为奇数" << endl;
return false;
}

int pause_4 = 0;
cout << "跳过直方图匹配" << endl;
cin >> pause_4;

if (!pause_4)
Histo_Match_Image(leftinfi, rightinfi);


int halfwindow = (windowsize) / 2;
Mat leftMat_Ori = imread(leftinfi.outputpath, CV_LOAD_IMAGE_ANYDEPTH);
Mat leftMat; leftMat_Ori.convertTo(leftMat, CV_64F);
Mat rightMat_Ori = imread(rightinfi.outputpath, CV_LOAD_IMAGE_ANYDEPTH);


Mat rightMat = rightMat_Ori;
int leftwidth_ = leftMat.cols;
int leftheight_ = leftMat.rows;
int rightwidth_ = rightMat.cols;
int rightheight_ = rightMat.rows;

Rect leftrange(0, 0, leftwidth_, leftheight_);
Rect rightrange(0, 0, rightwidth_, rightheight_);
vector<Point2d> _left_lsm;
vector<Point2d> _right_lsm;
vector<Point2d> _hfactor;
vector<Point2d> _leftmove;
vector<Point2d> _rightmove;
vector<double> _correvaluearray;
///对每个匹配点进行最小二乘精匹配
int icc = 0;
vector<Point2d> _leftnewkey;
vector<Point2d> _rightnewkey;
for (int ii = 0; ii < leftkey_.size(); ii++)
{
ResampleKernels* my_kernel = new ResampleKernels;
my_kernel->Init(Cubic_Convolution_6P);
int realwidth_ = windowsize - 2;
int realheight_ = windowsize - 2;
Rect arealeft(leftkey_[ii].x - halfwindow + 1, leftkey_[ii].y - halfwindow + 1, windowsize, windowsize);
Rect arearight(rightkey_[ii].x - halfwindow + 1, rightkey_[ii].y - halfwindow + 1, windowsize, windowsize);
if (!safeopencv::Inside(leftrange, arealeft)) continue; /// 此点在边界，暂时被忽略

///系数矩阵
double* coeffarray = new double[realwidth_ * realwidth_ * 8];
double* staticarray = new double[realwidth_ * realheight_]; // 观测值
double* tarray = new double[realwidth_ * realheight_ * 8];  // 常量矩阵

double CTC[64] = { 0 };
double CTL[8] = { 0 };

double pt[4] = { 0 };
double P0[4] = { leftkey_[ii].x, leftkey_[ii].y, rightkey_[ii].x, rightkey_[ii].y };

///左右片变形参数
double a0 = 0, a1 = 1, a2 = 0, b0 = 0, b1 = 0, b2 = 1;
///辐射畸变参数
double h0 = 0;
double h1 = 1;

///左影像数据、指针
Mat Left_Block_ori = leftMat(arealeft).clone();
Mat Left_Block; Left_Block_ori.convertTo(Left_Block, CV_64F);
void* pvoidleft = static_cast<void*>(Left_Block.data);
double* pimgleft = new double[windowsize * windowsize];// = static_cast<double*>(pvoidleft);

///左影像像平面坐标,真特么麻烦
vector<Point2f> Left_Block_Pixel;
safeopencv::Mat_Fill(Left_Block_Pixel, windowsize, Enum_Mat_Fill_Symmetric);

///填充左影像数据
for (int i = 0; i < windowsize * windowsize; i++)
{
Point2f leftresample_pos(Left_Block_Pixel[i].x + leftkey_[ii].x, Left_Block_Pixel[i].y + leftkey_[ii].y);
double tempvalue = safeopencv::Resample(rightMat, leftresample_pos, my_kernel);
pimgleft[i] = tempvalue;
}

///左、右窗口。  左窗口不动,右窗口初始值为粗匹配计算结果
///开始最小二乘参数的解算
double re[2] = { 0.95, -1 };
double th0, th1, ta0, ta1, ta2, tb1, tb2, tb0;
int ncounter = 0;
do
{
if (ncounter++ > 10)
break;
//判断右窗口是否仍在影像范围内
if (!safeopencv::Inside(rightrange, arearight)) {
cout << "右窗口溢出！！" << endl;
}
vector<Point2f> Right_Blcok_Pixel(windowsize * windowsize);
double a_affine_trans_data[6] = { a1, a2, a0,
b1, b2, b0 };
Mat affine_mat(2, 3, CV_32F, &a_affine_trans_data);
Point2f pixelstart;
pixelstart.x = rightkey_[ii].x; pixelstart.y = rightkey_[ii].y;
safeopencv::Pixel_Warp_Affine(Left_Block_Pixel, Right_Blcok_Pixel,pixelstart, a_affine_trans_data);


///重采样, 辐射畸变改正
Mat Right_Block(windowsize, windowsize, CV_64F);
int window = windowsize * windowsize;
void* pvoidright = static_cast<void*>(Right_Block.data);
double* pimgright = static_cast<double*>(pvoidright);



for (int q = 0; q < window; q++)
{
//cout << safeopencv::Resample(rightMat, Right_Blcok_Pixel[q], my_kernel) << endl;
double tempvalue = safeopencv::Resample(rightMat, Right_Blcok_Pixel[q], my_kernel);
pimgright[q] = h0 + h1 * tempvalue;
/*
int iix = Right_Blcok_Pixel[q].x;
int iiy = Right_Blcok_Pixel[q].y;
cout << safeopencv::Resample(rightMat, Right_Blcok_Pixel[q], my_kernel) << endl;
cout << rightMat.at<double>(iix, iiy) << endl;
cout << "111。。。。" << endl;
cout << iix << "\t" << iiy << endl;
int passs_1;
cin >> passs_1;

            }
            ///Right_Block_Pixel 为了方便重采样,此时是原影响坐标,需要转换
            ///计算相关系数
            re[1] = CalcCorrelation(pimgleft, pimgright, windowsize - 2, windowsize - 2);
            //cout << re[1] << endl;

            if (re[1]<re[0] && re[0]>0.95)
                break;

            re[0] = re[1];

            double a[8], aa[64], al[8];
            double x[8];
            cout.precision(4);
            Mat hei_test_left = leftMat(arealeft);
            Mat hei_test_right = rightMat(arearight);
            for (int j = 0; j < windowsize - 2; j++)
            {
                for (int i = 0; i < windowsize - 2; i++)
                {
                    int index_now = (j * (windowsize - 2) + i) * 8;
                    coeffarray[index_now + 0] = 1;
                    coeffarray[index_now + 1] = pimgleft[(j + 1) * windowsize + i + 1];
                    coeffarray[index_now + 2] = h1*(pimgright[(j + 1) * windowsize + i + 2] - pimgright[(j + 1) * windowsize + i])*0.5;
                    coeffarray[index_now + 3] = h1*(pimgright[(j + 1) * windowsize + i + 2] - pimgright[(j + 1) * windowsize + i])*0.5*Left_Block_Pixel[(j + 1) * windowsize + i + 1].x;//lx[i + 1][j + 1];
                    coeffarray[index_now + 4] = h1*(pimgright[(j + 1) * windowsize + i + 2] - pimgright[(j + 1) * windowsize + i])*0.5*Left_Block_Pixel[(j + 1) * windowsize + i + 1].y;//ly[i + 1][j + 1];
                    coeffarray[index_now + 5] = h1*(pimgright[(j + 2) * windowsize + i + 1] - pimgright[(j + 0) * windowsize + i + 1])*0.5;
                    coeffarray[index_now + 6] = h1*(pimgright[(j + 2) * windowsize + i + 1] - pimgright[(j + 0) * windowsize + i + 1])*0.5*Left_Block_Pixel[(j + 1) * windowsize + i + 1].x;
                    coeffarray[index_now + 7] = h1*(pimgright[(j + 2) * windowsize + i + 1] - pimgright[(j + 0) * windowsize + i + 1])*0.5*Left_Block_Pixel[(j + 1) * windowsize + i + 1].y;

                    staticarray[j*(windowsize - 2) + i] = pimgleft[(j + 1) * windowsize + i + 1] - pimgright[(j + 1) * windowsize + i + 1];


                
                    cout << "calced data:" << endl;
                    cout << pimgleft[(j + 1) * windowsize + i + 1] << "\t" << pimgright[(j + 1) * windowsize + i + 1] << endl;
                    cout << leftMat.at<double>(leftkey_[ii].y - halfwindow + 1 + i + 1, leftkey_[ii].x - halfwindow + j + 1) << "\t" << rightMat.at<double>(rightkey_[ii].y - halfwindow + 1 + i + 1, rightkey_[ii].x - halfwindow + j + 1) << "\t" << endl;
                    
                }

            }


            ///使用逐点法化
            helper.transpose(coeffarray, tarray, realwidth_ * realheight_, 8);
            helper.mult(tarray, coeffarray, CTC, 8, realwidth_*realheight_, 8);
            helper.mult(tarray, staticarray, CTL, 8, realwidth_*realheight_, 1);



            //CTL为8X1
            if (!helper.invers_matrix(CTC, 8))		 					//求逆
            {
                cout << "不可逆" << endl;
                ;
            }

            double V[8];


            helper.mult(CTC, CTL, V, 8, 8, 1);
            memcpy(x, al, sizeof(double)* 8);
            cout.precision(20);


            //V对应dh0 dh1  da0-da2 db0-db2
            th0 = h0; th1 = h1; ta0 = a0; ta1 = a1; ta2 = a2; tb1 = b1; tb2 = b2; tb0 = b0;
            h0 = th0 + V[0] + th0*V[1];
            h1 = th1 + th1*V[1];
            a0 = ta0 + V[2] + ta0*V[3] + tb0*V[4];
            a1 = ta1 + ta1*V[3] + tb1*V[4];
            a2 = ta2 + ta2*V[3] + tb2*V[4];
            b0 = tb0 + V[5] + ta0*V[6] + tb0*V[7];
            b1 = tb1 + ta1*V[6] + tb1*V[7];
            b2 = tb2 + ta2*V[6] + tb2*V[7];


        } while (true);


        double sum_x = 0, sum_y = 0, sum_gx = 0, sum_gy = 0;
        for (int j = 0; j < realheight_; j++)
        {
            for (int i = 0; i < realwidth_; i++)
            {
                sum_x += Left_Block_Pixel[i + 1 + (j + 1) *windowsize].x * (pimgleft[i + 2 + (j + 1) * windowsize] - pimgleft[i + (j + 1) * windowsize])*(pimgleft[i + 2 + (j + 1) * windowsize] - pimgleft[i + (j + 1) * windowsize])*0.25;
                sum_y += Left_Block_Pixel[i + 1 + (j + 1) *windowsize].y * (pimgleft[i + 1 + (j + 2) * windowsize] - pimgleft[i + 1 + j * windowsize])*(pimgleft[i + 1 + (j + 2) * windowsize] - pimgleft[i + 1 + j * windowsize])*0.25;
                sum_gx += (pimgleft[i + 2 + (j + 1) * windowsize] - pimgleft[i + (j + 1) * windowsize])*(pimgleft[i + 2 + (j + 1) * windowsize] - pimgleft[i + (j + 1) * windowsize])*0.25;
                sum_gy += (pimgleft[i + 1 + (j + 2) * windowsize] - pimgleft[i + 1 + j * windowsize])*(pimgleft[i + 1 + (j + 2) * windowsize] - pimgleft[i + 1 + j * windowsize])*0.25;
            }
        }

        double plusleft1 = pt[0] = sum_x / sum_gx;
        double plusleft2 = pt[1] = sum_y / sum_gy;
        double plusright1 = pt[2] = ta0 + ta1*pt[0] + ta2*pt[1];
        double plusright2 = pt[3] = tb0 + tb1*pt[0] + tb2*pt[1];


        //if (re[1] > 0.5)
        //cout << "pt:" << pt[0] << "\t" << pt[1] << "\t" << pt[2] << "\t" << pt[3] << endl;
        pt[0] += P0[0];
        pt[1] += P0[1];
        pt[2] += P0[2];
        pt[3] += P0[3];
        pt[4] = re[0];



        //if (re[1] > 0.5)

        //cout << "matches point:" << pt[0] << "\t" << pt[1] << "\t" << pt[2] << "\t" << pt[3] << endl;
        if (re[1] > 0.9)
        {
            _left_lsm.push_back(Point2d(pt[0], pt[1]));
            _right_lsm.push_back(Point2d(pt[2], pt[3]));
            _hfactor.push_back(Point2d(th0, th1));
            _leftmove.push_back(Point2d(plusleft1, plusleft2));
            _rightmove.push_back(Point2d(plusright1, plusright2));
            _correvaluearray.push_back(re[1]);

            double rateleft = 0.9 / sqrt(pow(plusleft1, 2.0) + pow(plusleft2, 2.0));
            _leftnewkey.push_back(Point2d(pt[0] + rateleft * plusleft1, pt[1] + rateleft * plusleft2));
            double rateright = 0.9 / sqrt(pow(plusright1, 2.0) + pow(plusright2, 2.0));
            _rightnewkey.push_back(Point2d(pt[2] + rateright * plusright1, pt[3] + rateright * plusright2));
        }


    }
*/

bool ImgFacotryAPI::Form_Template()
{
    Template_PreMatch("c:\\Temp\\match\\heimatch.xls");

    int startindex_x = 20;
    int startindex_y = 20;

    int endindex_x = 459;
    int endindex_y = 10680;

    _leftkey_.clear();
    _rightkey_.clear();

    for (int j = startindex_y; j < endindex_y; j += 5)
    {
        for (int i = startindex_x; i < endindex_x; i += 5)
        {
            double xvalue = barycentriccalc(x2xinfo, static_cast<double>(i));
            double yvalue = barycentriccalc(y2yinfo, static_cast<double>(j));

            if (xvalue < 20 || xvalue > 459 || yvalue < 20 || yvalue > 10680)
                continue;

            else
            {
                _leftkey_.push_back(Point2d(static_cast<double>(i), static_cast<double>(j)));
                _rightkey_.push_back(Point2d(xvalue,yvalue));
            }

        }
    }

    return true;
}

bool ImgFacotryAPI::Match_Template(MatchToolType leftinfi, MatchToolType rightinfi, int windowsize)
{
    lsm_multi_times = false;
    string img_left_path = leftinfi.filepath;
    string img_right_path = leftinfi.filepath;

    Mat left_mat = imread(img_left_path, IMREAD_ANYDEPTH);
    Mat right_mat = imread(img_right_path, IMREAD_ANYDEPTH);

    if (windowsize % 2 == 0)
    {
        cout << "请将模板匹配的窗口大小设置为奇数" << endl;
        return false;
    }

    int bigwindow = (windowsize)* 2 - 1;

    vector<Point2d> left_result;
    vector<Point2d> right_result;
    for (int i = 0; i < _leftkey_.size(); i++)
    {
        Mat template_content = left_mat(Rect(_leftkey_[i].x - windowsize / 2, _leftkey_[i].y - windowsize / 2, windowsize, windowsize)).clone();
        Mat range_content = left_mat(Rect(_rightkey_[i].x - bigwindow / 2, _rightkey_[i].y - bigwindow / 2, bigwindow, bigwindow)).clone();
    
        //模版匹配
        /// 创建输出结果的矩阵
        int result_cols = range_content.cols - template_content.cols + 1;
        int result_rows = range_content.rows - template_content.rows + 1;
        Mat result;
        result.create(result_cols, result_rows, CV_32FC1);

        /// 进行匹配和标准化
        matchTemplate(range_content, template_content, result, CV_TM_CCOEFF_NORMED);
        normalize(result, result, 0, 1, NORM_MINMAX, -1, Mat());

        double minVal; double maxVal; Point minLoc; Point maxLoc;
        Point matchLoc;

        minMaxLoc(result, &minVal, &maxVal, &minLoc, &maxLoc, Mat());


        left_result.push_back(_leftkey_[i]);
        right_result.push_back(Point2d(_rightkey_[i].x + maxLoc.x - windowsize/2, _rightkey_[i].y + maxLoc.y - windowsize/2 ));
    }

    _leftkey_.clear();
    _rightkey_.clear();

    for (int i = 0; i < left_result.size(); i++)
    {
        _leftkey_.push_back(left_result[i]);
        _rightkey_.push_back(right_result[i]);
    }

    cout << "模版匹配得到的初值" << endl;
    cout << _leftkey_.size() << endl;

    cout << "进行模板匹配过程中的最小二乘精匹配？" << endl;
    int match_mode = 0;
    cin >> match_mode;
    if (match_mode == 1)
    {
        LSM_Match(leftinfi, rightinfi, 19);
    }
    return false;
}


