#include "stdafx.h"
#include <new>
#include "ImgFacotryAPI.h"
//#define SKIP_PROCESS_IMAGEOUT
#define CONDUCT_OUTERIOR_PARA_CALIBRATION_PROCESS_
//#define _ImageOutPutCheckingStart_Mode_Test

bool VirtualMachine::GetAuxFile(string _auxfile, CEnumMachineCompent _COMPonetType)
{
    ifstream ifs;
    ifs.open(_auxfile.c_str(), ios::in);
    std::stringstream buffer;
    buffer << ifs.rdbuf();
    std::string contents(buffer.str());

    string realcontentes = safestring::SingleLine(contents);
    map<CEnumMachineCompent, string>::iterator itr;
    itr = _compent_lookuptable.find(_COMPonetType);
    char tempreader[1024];
    sprintf_s(tempreader, 1024, "<%s>(([\\s\\S]*?))</%s>", itr->second.c_str(), itr->second.c_str());
    //(([\s\S]*?))
    vector<string> orbitapistring = safestring::GetRegexResult(realcontentes, static_cast<string>(tempreader));
    for (int i = 0; i < orbitapistring.size(); i++)
    {
        //cout << (static_cast<double>(i)) / (orbitapistring.size())*100 << "\r";
        if (_COMPonetType == MachineComponent_asAttitude)
        {
            AttitudeAPI tempapi(orbitapistring[i], asXml);
            _attitudeapi.push_back(tempapi);
        }
        else if (_COMPonetType == MachineComponent_asLienTime)
        {
            LineTimeAPI tempapi(orbitapistring[i], asXml);
            _timeapi.push_back(tempapi);
        }
        else if (_COMPonetType == MachineComponent_asOrbit)
        {
            OrbitAPI tempapi(orbitapistring[i], asXml);
            _orbitApi.push_back(tempapi);
        }
    }

    ifs.close();
    return true;
}

UTCAPI VirtualMachine::GetTime(CEnumMachineCompent _apiname, int _indexid)
{
    string _timestringofit;
    if (_apiname == MachineComponent_asAttitude)
        _timestringofit = _attitudeapi[_indexid].IGetTime();

    else if (_apiname == MachineComponent_asLienTime)
        _timestringofit = _timeapi[_indexid].IGetTime();
    else if (_apiname == MachineComponent_asOrbit)
        _timestringofit = _orbitApi[_indexid].IGetTime();

    UTCAPI utcapi(_timestringofit, Enum_TimeStandard_Standard);
    return utcapi;
}

void VirtualMachine::SetPropertys()
{
    char tempreader[1024];
    cout << "请输入轨道根节点的命名:" << endl;
    cin >> tempreader;
    _compent_lookuptable[MachineComponent_asOrbit] = tempreader;

    cout << "请输入姿态根节点的命名：" << endl;
    cin >> tempreader;
    _compent_lookuptable[MachineComponent_asAttitude] = tempreader;

    cout << "请输入时间根节点的命名:" << endl;
    cin >> tempreader;
    _compent_lookuptable[MachineComponent_asLienTime] = tempreader;
}

void VirtualMachine::SetLookUpTable()
{
    _compent_lookuptable[MachineComponent_asOrbit] = "GpsData";
    _compent_lookuptable[MachineComponent_asAttitude] = "AttData";
    _compent_lookuptable[MachineComponent_asLienTime] = "LineTime";
}

VirtualMachine::VirtualMachine(string auxxml, bool withimage) : base_point(SetPoint())
{
    SetLookUpTable();
    if (auxxml == "UnitTest")
    {
        SetEarthModel();
        MetaClass::MetaClass("auxconfig");
        return;
    }
    GetAuxFile(auxxml, MachineComponent_asLienTime);
    GetAuxFile(auxxml, MachineComponent_asAttitude);
    GetAuxFile(auxxml, MachineComponent_asOrbit);
}

bool VirtualMachine::GetJSonFile(string _filepath, CEnumMachineCompent _componenttype)
{
    if (_componenttype == MachineComponent_asAttitude)
        _filepath += "\\" + string("AttitudeInfo.json");
    else if (_componenttype == MachineComponent_asLienTime)
        _filepath += "\\" + string("LineTimeinfo.json");
    else if (_componenttype == MachineComponent_asOrbit)
        _filepath += "\\" + string("OrbitInfo.json");
    ifstream ifs(_filepath.c_str());
    char reader[2048];
    do
    {
        ifs.getline(reader, 1024);
        if (_componenttype == MachineComponent_asAttitude)
        {
            AttitudeAPI attitudeapi(reader, asJSon);
            _attitudeapi.push_back(attitudeapi);
        }
        else if (_componenttype == MachineComponent_asLienTime)
        {
            LineTimeAPI linetimeapi(reader, asJSon);
            _timeapi.push_back(linetimeapi);
        }
        else if (_componenttype == MachineComponent_asOrbit)
        {
            OrbitAPI orbitapi(reader, asJSon);
            _orbitApi.push_back(orbitapi);
        }
    } while (ifs.peek() != EOF);

    return true;
}

VirtualMachine::VirtualMachine(string dirpath, CEnumDomType _jsonobject) : base_point(SetPoint())
{
    SetLookUpTable();
    MetaClass::MetaClass("auxconfig");
    GetJSonFile("d:\\testfactory\\jsonhost", MachineComponent_asAttitude);
    GetJSonFile("d:\\testfactory\\jsonhost", MachineComponent_asLienTime);
    GetJSonFile("d:\\testfactory\\jsonhost", MachineComponent_asOrbit);

}

int VirtualMachine::GetSize(CEnumMachineCompent enu)
{
    if (enu == MachineComponent_asAttitude)
    {
        return _attitudeapi.size();
    }
    else if (enu == MachineComponent_asLienTime)
    {
        return _timeapi.size();
    }

    else if (enu == MachineComponent_asOrbit)
    {
        return _orbitApi.size();
    }
}

string VirtualMachine::IGetPosition(double _utctime)
{
    OrbitAPI hei = _orbitApi[0];
    m_base.LagrangianInterpolation(_orbitApi, _utctime, _orbitApi.size(), _orbitApi[0]);
    return _orbitApi[0].Serialize();
}

AttitudeAPI VirtualMachine::GetAttitude(double _utctime)
{
    AttitudeAPI hei = _attitudeapi[0];
    m_base.QuatInterpolation(_attitudeapi, _utctime, _attitudeapi.size(), hei);
    return hei;
}

/*
int lineid;
int actLine;
double utctime;
double avgtime;
*/
string VirtualMachine::OutPutTime()
{
    ofstream ofs;
    string filepath = "TimeCheck" + safestring::WholeTimeString() + ".xls";
    int sizeoftime = _simpleTimes.size();
    ofs.open(filepath.c_str(), ios::out);
    ofs << "lineid" << "\t" << "actLine" << "\t" << "UTCTime" << "\t" << "AvgTime" << endl;
    ofs.precision(20);
    for (int i = 0; i < sizeoftime; i++)
    {
        ofs << _simpleTimes[i].lineid << "\t" << _simpleTimes[i].actLine << "\t" << _simpleTimes[i].utctime
            << "\t" << _simpleTimes[i].avgtime << endl;
    }

    ofs.close();

    return filepath;
}

bool VirtualMachine::BuildRpcModel(int index)
{
    int imagesizeofit[4] = { 0, 0, _singleframeWidth, _singleframeHeight };
    _rpcInfo[index - 41].CreateRPCModel(this, true, imagesizeofit, index, 3, 3, 5);

    return true;
}

void  VirtualMachine::FromXY2LONLATRPC(const double& imgx, const double& imgy, const double& height, double& lon, double&lat, const int& index)
{
    _rpcInfo[index - 41].FromImageToGround(imgx, imgy, height, lon, lat);
}

double VirtualMachine::GetAlongAngle(double yvalue,bool bxarray)
{
    double getvalue = 0;
    if (!bxarray)
    {
        vector<TimeCalibStruct>::iterator itr = p_TimeCalib->begin();

        for (itr = p_TimeCalib->begin(); itr != p_TimeCalib->end(); itr++)
        {
            if (yvalue > itr->_startvalue && yvalue < itr->_endvalue)
            {
                double intervalue = yvalue;
                getvalue = barycentriccalc(itr->_interkernel, intervalue);
                break;
            }

            else if ((itr + 1) != p_TimeCalib->end() && yvalue > itr->_endvalue && yvalue < (itr + 1)->_startvalue)
            {
                double calleft = barycentriccalc(itr->_interkernel, itr->_endvalue);
                double calright = barycentriccalc((itr + 1)->_interkernel, (itr + 1)->_startvalue);

                double distance = (itr + 1)->_startvalue - itr->_endvalue;
                getvalue = calleft * ((itr + 1)->_startvalue - yvalue) / distance + calright * (yvalue - itr->_endvalue) / distance;
                break;
            }
        }

    }
    else
    {
        vector<TimeCalibStruct>::iterator itr = p_TimeCalib_X->begin();

        //cout <<"interpolant value is :"<< yvalue << endl;
        for (itr = p_TimeCalib_X->begin(); itr != p_TimeCalib_X->end(); itr++)
        {
            if (yvalue > itr->_startvalue && yvalue < itr->_endvalue)
            {
                double intervalue = yvalue;
                getvalue = barycentriccalc(itr->_interkernel, intervalue);
                break;
            }

            else if ((itr+1) != p_TimeCalib_X->end() && yvalue > itr->_endvalue && yvalue < (itr + 1)->_startvalue)
            {
                double calleft = barycentriccalc(itr->_interkernel, itr->_endvalue);
                double calright = barycentriccalc((itr + 1)->_interkernel, (itr + 1)->_startvalue);

                double distance = (itr + 1)->_startvalue - itr->_endvalue;
                getvalue = calleft * ((itr + 1)->_startvalue - yvalue) / distance + calright * (yvalue - itr->_endvalue) / distance;
                break;
            }
        }
    }
    //cout << getvalue << endl;
    return getvalue;
}

VirtualMachine::VirtualMachine(string fredconfig, string fredpath) :base_point(SetPoint())
{
	rpc_container.reserve(100);
	rpc_container.resize(100);
    p_TimeCalib = nullptr;
    p_TimeCalib_X = nullptr;
    _rpcInfo.resize(100);
    GDALAllRegister();
    EnvelopeAPI envAPI(fredconfig, fredpath);

    if (true/*envAPI.IFindStart() == Process_Success*/)
    {
        envAPI.IScan();
        //envAPI.ICheckContent();
        unsigned int imgframewidth = envAPI.IWidth();  _singleframeWidth = imgframewidth;
        unsigned int imgframeheight = envAPI.IHeight(); _singleframeHeight = imgframeheight;

        ///BASE_POINT has not been initialized yet
        base_point = SetPoint(imgframeheight, imgframewidth);
        unsigned int bytecount = envAPI.IByteCount();

        unsigned int framecount = envAPI.IFrameNum();
        unsigned int auxcount;
        cout << "Fred文件装载完毕，开始分离数据" << endl;


        ///分离出辅助数据
        cout << "开始分离辅助数据..." << endl;
        cout << "准确获取辅助数据注册文件，如果缺失，会有自动注册功能" << endl;
        //此处把配置文件的路径设置固定了 以后要进行修改
        MetaClass::MetaClass("auxconfig");
        AuxFactory auxfct("auxBinaryconfig");
        auxfct.ClearFrame();
        auxcount = envAPI.IAuxNum();

        int datainfo = framecount / auxcount; int framefindstart = 40;
        bool findauxstart = false;
        for (int i = 0; i < auxcount; i++)
        {
            AuxFactory myaxu = envAPI.IFormAux();
            //envAPI.ICheckContent();
            //检查数据正确性 ，然后对辅助数据更新 
            if (myaxu.Attitude().IGetTime() == 0)
                continue;
            else
            {
                if (!findauxstart)
                {
                    //framefindstart = i * datainfo;
                    findauxstart = true;

                }
                int framelimit = (11265.0) * 0.5 *1.2;
                //排除意外情况

                if (_myFactory.size() >= 1)
                {
                    if (fabs(static_cast<double>(myaxu.GetCol()) - _myFactory[_myFactory.size() - 1].GetCol()) <= 1)
                    {
                        cout << "异常，重现列数重复辅助数据，序列号：" << i << "和" << i + 1 << endl;
                        IAddLog("异常，重现列数重复辅助数据，序列号：%d和%d", i, i + 1);
                    }
                    if ((static_cast<double>(myaxu.GetCol()) - _myFactory[_myFactory.size() - 1].GetCol()) > 0)
                    {

                    }
                }

                myaxu.AddFrame();
                myaxu.FrameNum(myaxu.framecount);
                _myFactory.push_back(myaxu);
            }

        }

        cout << "有效的辅助数据帧数" << _myFactory.size() << endl;
        cout << "有效的辅助数据对应的景数（未顾及跳变)" << _myFactory[0].framecount << endl;
        IAddLog("有效的辅助数据帧数:%d;有效的辅助数据对应的景数:%d（未顾及跳变).", _myFactory.size(), _myFactory[0].framecount);
        cout << "跳变数据已失效" << endl;
        cout << "开始建立成像框架" << endl; IAddLog("开始建立成像框架");
        m_pImgModel = new ImageModel(FormImageModel());
        //自检
        if (!m_pImgModel->SelfCheck())
        {
            cout << "组建的模型自检不通过，请检查日志" << endl;
            return;
        }
        cout << "成像框架建立完毕" << endl; IAddLog("成像框架建立完毕");
        SetEarthModel();




        ///整理行时 
        ///构建，优化行时
        cout << "优化行时数据中" << endl;
        IAddLog("优化行时数据中");
        vector<int> heitemp;
        heitemp.push_back(_singleframeHeight);
        SetSimpleTime(heitemp);
        //将行时以EXCEL规定的格式输出
        //string filepathftime = OutPutTime();
        //cout << "简化的行时数据被输出在" << filepathftime.c_str() << "中以供检查测试" << endl;
        //cout << "行时数据优化完毕..." << endl;
        IAddLog("行时数据优化完毕...");


        base_point.Serialize("testbasepoint.json");
        /*
        基本 定位 测试
        */
        double testx, testy;
        testx = 5070; testy = 0;
        double height = 100;
        double lontest, latest;
        for (int i = 0; i < 47; i++)
        {
            for (int j = 0; j < (10786 - 10); j += 85)
                ;// FromXY2LonLatTest(i * 10, j, height, lontest, latest, 62);
        }

        int container_i[10] = { 62, 65, 66, 67, 72, 74, 75, 76, 78, 79 };
        for (int i = 0; i < 10; i++)
        {
            int index_i = container_i[i];
            FromXY2LonLatTest(0, 0, 100, lontest, latest, index_i);
            double oldtestlon, oldtestlat;
            oldtestlon = lontest;
            oldtestlat = latest;
            cout << index_i - 61 << "(0,0)->(" << lontest << "," << latest << ")" <<  endl;
            FromXY2LonLatTest(0, 10785, 100, lontest, latest, index_i);
            cout << index_i - 61 << "(0,10785)->(" << lontest << "," << latest << ")" << "distatnce:("<<lontest-oldtestlon<<","<<latest - oldtestlat<<")" <<endl;
        }





        //ReadSwingError();

#ifdef SKIP_PROCESS_IMAGEOUT
        cout << "程序即将结束" << endl;
        return;
#endif

#ifdef CONDUCT_OUTERIOR_PARA_CALIBRATION_PROCESS_

        ifstream readcontrolpoint;
        readcontrolpoint.open("C:\\Temp\\CONTROLPOINTS.json", ios::in);
        if (!readcontrolpoint.is_open())
        {
            cout << "控制点文件不存在,控制点外畸变参数校正将不被执行,详情请见日志" << endl;
            IAddLog("控制点文件%s不存在，请检查配置，相应的控制点外畸变校正工作将不被畸形，畸变参数将不会与之前发生变化", "C:\\Temp\\CONTROLPOINTS.json");
        }
        else
        {
            vector<ImgPoint> controlist;
            do
            {
                char readcontrolLine[1024];
                readcontrolpoint.getline(readcontrolLine, 1024);
                //需要行时 
                ImgPoint newpoints(static_cast<string>(readcontrolLine), asJSon);
                int _lineindex = currentImageID * _singleframeHeight + _singleframeHeight - newpoints.ImgY();
                double time = _simpleTimes[_lineindex].utctime;  //时间的设置 
                newpoints.utcapi = UTCAPI(time, Enum_TimeStandard_ZY03);
                newpoints = newpoints + base_point;
                double xvalue, yvalue, zvalue;
                double heilon, heilat;
                FromXY2LonLatTest(newpoints.ImgX(), newpoints.ImgY(), 100, heilon, heilat, 62);
                pEarth->Geo2Rect(-heilon* CV_PI / 180.0, heilat * CV_PI / 180.0, 100, xvalue, yvalue, zvalue);
                double xposvalue, yposvalue, zposvalue;
                pEarth->Geo2Rect(-newpoints.Lon() * CV_PI / 180.0, newpoints.Lat() * CV_PI / 180.0, 100, xposvalue, yposvalue, zposvalue);

                double test1, test2, test3;
                newpoints.OriX(xvalue);
                newpoints.OriY(yvalue);
                newpoints.OriZ(zvalue);
                newpoints.X(xposvalue); newpoints.Y(yposvalue); newpoints.Z(zposvalue);
                newpoints.Serialize("d:\\testPos\\0410\\");
                double heisqrt = 0;
                heisqrt += pow(newpoints.X() - newpoints.OriX(), 2.0);
                heisqrt += pow(newpoints.Y() - newpoints.OriY(), 2.0);
                heisqrt += pow(newpoints.Z() - newpoints.OriZ(), 2.0);
                heisqrt = sqrt(heisqrt);
                cout << heisqrt << "\t" << heilon << "\t" << heilat << "\t" << newpoints.Lon() << "\t" << newpoints.Lat() << endl;
                if (fabs(newpoints.Lat() - 0) < 0.001)
                    continue;
                controlist.push_back(newpoints);
            } while (readcontrolpoint.peek() != EOF);

            readcontrolpoint.close();
            CalibratedModel calibration(m_pImgModel);
            calibration.Calibration(controlist,base_point);

            ifstream ifs_swing("c:\\Temp\\match\\heimatch_swing.xls");
            ofstream ofs_swing("c:\\Temp\\match\\heimatch_swing_error_calc.xls");
            int ncounter_swing = 0;
            ofs_swing.precision(20);
            double add_number = 0;
            do
            {
                char read_swing[1024];
                ifs_swing.getline(read_swing, 1024);
                double a_swing[6];
                sscanf_s(read_swing, "%lf %lf %lf %lf %lf", &a_swing[0], &a_swing[1], &a_swing[2], &a_swing[3], &a_swing[4]);
                add_number += a_swing[4];

                ncounter_swing++;
                a_swing[5] = add_number - ncounter_swing;

                char printf_swing[1024];
                sprintf_s(printf_swing, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf", a_swing[0], a_swing[1], a_swing[2], a_swing[3], a_swing[4], a_swing[5]);
                ofs_swing << printf_swing << endl;
            } while (ifs_swing.peek() != EOF);
            ifs_swing.close();
            ofs_swing.close();

            cout << "进行外检校?" << endl;
            int calibmode = 0;
            cin >> calibmode;
            if (calibmode == 1)
            {
                CalibratedModel calibrate(m_pImgModel);
                vector<ImgPoint> allve;
                GetControlPoints("C:\\Temp\\match\\heimatch.xls",allve);
                cout << endl;
                cout << "控制点读取完毕" << endl;
                calibrate.Calibration(allve,base_point);
            }



            cout << "进行时间标定和摆扫角标定？" << endl;
            int calibmode_2 = 0;
            cin >> calibmode_2;
            if (calibmode_2 == 1)
            {
                int cuid = 62;
                double crabangle = CalcCrabAngle(cuid);
                double Vx = Vel(cuid);
                double kappaOfCamera = 1.4 * CV_PI / 180.0;
                double rangestep = GetRangeForOneStep(62);
                CalibratedModel calibrate(m_pImgModel);
                vector<ImgPoint> allve;
                GetControlPoints("C:\\Temp\\match\\heimatch.xls", allve);
                cout << endl;

                vector<SwingCablib> heiswing;

                double TimeSegMent = perT * _singleframeHeight;
                calibrate.CalibrateSwing(allve, base_point, _gofuckyourselfSwingError,endswingvalue);
                double endxtimevalue;
                calibrate.CalibrateTime_(allve, base_point, heiswing, endxtimevalue);
            }
            
            cout << "进行摆扫度非线性标定?" << endl;
            int calibmode_swing_linear = 0;
            cin >> calibmode_swing_linear;
            if (calibmode_swing_linear == 1)
            {
                CalibratedModel calibrate(m_pImgModel);
                vector<ImgPoint> allve;
                GetControlPoints("C:\\Temp\\match\\heimatch.xls", allve,true);
                calibrate.Calibration_Swing(allve, base_point, 20, 10670, 5);
            }

            cout << "进行指向角标定?" << endl;
            int calibmode_3 = 0;
            cin >> calibmode_3;
            if (calibmode_3 == 1)
            {
                cout << "进行匹配曲线拟合?" << endl;
                int curvematch_mode = 0;
                cin >> curvematch_mode;
                
                if (curvematch_mode == 1)
                {
                    ImgFacotryAPI imgapitemp(480, 10786, 12, "hei");
                    imgapitemp.Form_Template();
                    MatchToolType left_type;
                    MatchToolType right_type;
                    left_type.filepath = left_type.outputpath = "C:\\Temp\\match\\left_.tif";
                    right_type.filepath = right_type.outputpath = "C:\\Temp\\match\\right_PreMatch.tiff8bit.tif";
                    
                    left_type.outputpath = "D:\\OpenCVWithCuda\\RPCModel\\ConsoleApplication1\\ConsoleApplication1\\Whole.tif25_.tif";
                    right_type.outputpath = "C:\\Temp\\match\\right_PreMatch.tiffhei.tif";// "D:\\IRS\\JB11-1_IRS_000116063_003_002_L0\\JB11-1_IRS_000116063_003_004_002_L3C_TOABT.tiffPreMatch.tiff";
                    imgapitemp.Match_Template(left_type, right_type, 19);

                    cout << "匹配曲线拟合完毕: 是否继续?" << endl;
                    cin >> curvematch_mode;
                }

                ///察看指向角标定文件是否存在
                ifstream ifs_lookangle("C:\\Temp\\match\\lookangle_find_error.txt");
                ifstream ifs_lookangle_distance("C:\\Temp\\match\\look_error_test.txt");
                CalibratedModel calibrate(m_pImgModel);
                vector<ImgPoint> allve;
                GetControlPoints("C:\\Temp\\match\\heimatch.xls", allve);
                cout << endl;
                cout << "控制点读取完毕" << endl;
                calibrate.Calibration(allve, base_point);
                if (ifs_lookangle.is_open())
                {
                    cout << "指向角标定文件存在! 读取中..." << endl;
                    do
                    {
                        char readlook[1024];
                        char readerror[1024];
                        ifs_lookangle.getline(readlook, 1024);
                        ifs_lookangle_distance.getline(readerror, 1024);
                        int indexlook; double x_look, y_look;
                        int indexlook_test; double x_error, y_error;
                        sscanf_s(readlook, "%d %lf %lf", &indexlook, &x_look, &y_look);
                        sscanf_s(readerror, "%d %lf %lf", &indexlook_test, &x_error, &y_error);
                        if (indexlook != indexlook_test)
                        {
                            cout << "错误,指向角标定文件不完整!" << endl;
                            cout << "(" << indexlook << "," << indexlook_test << ")" << endl;
                            cout << "是否中止程序" << endl;
                            int pause_exit = 0;
                            cin >> pause_exit;
                            if (pause_exit == 1)
                                exit(0);
                        }
                        if (fabs(y_error) > 4) continue;
                        //if (indexlook == 84) continue;
                        
                        base_point.SetAngle(indexlook, x_look, y_look);
                    } while (ifs_lookangle.peek() != EOF);

                    ifs_lookangle.close();
                    cout << "指向角文件读取完毕" << endl;
                    CalibratedModel calibrate1(m_pImgModel);

                    vector<ImgPoint> allve1;
                    GetControlPoints("C:\\Temp\\match\\heimatch.xls", allve1);
                    cout << endl;
                    cout << "控制点读取完毕" << endl;
                    calibrate1.Calibration(allve1, base_point);
                }

                else
                {

                    int cuid = 62;
                    double crabangle = CalcCrabAngle(cuid);
                    double Vx = Vel(cuid);
                    double kappaOfCamera = 1.4 * CV_PI / 180.0;
                    double rangestep = GetRangeForOneStep(62);
                    CalibratedModel calibrate(m_pImgModel);
                    vector<ImgPoint> allve;
                    GetCompensate();
                    GetControlPoints("C:\\Temp\\match\\heimatch.xls", allve);
                    cout << endl;

                    vector<SwingCablib> heiswing;
                    calibrate.CalibrateLookAngle(allve, base_point);
                }
            }
        }
#endif 
        /*
        BuildRpcModel(62);
        double maxerrorx = -1000,maxerrory = -1000;
        double rmserrorx =  0, rmserrory = 0;
        int k = 0;
        for (int j = 0; j < _singleframeHeight; j+=5)
        {
            for (int i = 0; i < _singleframeWidth; i+=5,k++)
            {
                double lon1, lat1;
                double lon2, lat2;
                FromXY2LonLatTest(i, j, 100, lon1, lat1, 62);
                FromXY2LONLATRPC(i, j, 100, lon2, lat2, 62);
                maxerrorx = (fabs(lon1 - lon2) > maxerrorx) ? fabs(lon1 - lon2) : maxerrorx;
                maxerrory = (fabs(lat1 - lat2) > maxerrory) ? fabs(lat1 - lat2) : maxerrory;
                rmserrorx += pow(fabs(lon1 - lon2), 2.0);
                rmserrory += pow(fabs(lat1 - lat2), 2.0);
            }
        }
        rmserrory = sqrt(rmserrory /k -1);
        rmserrorx = sqrt(rmserrorx / k - 1);
        std::cout.precision(20);
        std::cout << "Maxerror between Rigurous Model and RPC" << "( " << maxerrorx << "," << maxerrory <<")" <<endl;
        std::cout << "RMS ERROR BETWEEN RIGUROUS MODEL AND RPC" << "( " << rmserrorx << "," << rmserrory << ")" << endl;
        return;
        */
        double lonres1, latres1;
        CalibratedModel calibration(m_pImgModel);
        

        ///进行计算闭合性检查 
		/*
        double circle_check_lonlat[2];
		double error_avg = 0;
		vector<double> record_error;
		int t_count = 0;
		for (int t_j = 0; t_j < _singleframeHeight - 1; t_j+=5)
		{
			cout << t_j << "\\" << _singleframeHeight - 1 << "\r";
			for (int t_i = 0; t_i < _singleframeWidth - 1; t_i += 5, t_count++)
			{
				double d_i = t_i; double d_j = t_j;
				//cout << d_i << "," << d_j << endl;
				FromXY2LonLatTest(d_i, d_j, 100, circle_check_lonlat[0], circle_check_lonlat[1], 66);
				double circle_check_xy[2];
				FromlatlonH2xy(circle_check_lonlat[0], circle_check_lonlat[1], 100, circle_check_xy[0], circle_check_xy[1], 66);
				double  error_cals = sqrt(pow(d_i - circle_check_xy[0], 2.0) + pow(d_j - circle_check_xy[1], 2.0));
				record_error.push_back(error_cals);
				if (fabs(error_cals) >1)
				{
					cout << "残差出错！" << endl;
					cout << d_i << "," << d_j << endl;
					int pause_mode_continue;
					cin >> pause_mode_continue;
				}
				error_avg += error_cals;
			}
		}

		error_avg /= t_count;
		double error_rms = 0;
		for (int t_i = 0; t_i < t_count; t_i++)
		{
			error_rms += pow(record_error[t_i] - error_avg, 2.0) / (t_count -1);
		}
  
		error_rms = sqrt(error_rms);


		cout << "正反算平均误差：" << error_avg;
		cout << "正反算残差中误差:" << error_rms;
        cout << "是否继续" << endl;
        int pause_mode_continue;
        cin >> pause_mode_continue;
		*/

		double lon_lat_mosaic[2];
		double xy_mosaic[2];


		
		cout << "生成rpc?：" << endl;
		int rpc_gener = 0;

		ImgFacotryAPI imgapi1(imgframewidth, imgframeheight, bytecount, "hei");
		cin >> rpc_gener;
		if (rpc_gener == 1)
		{ 
			Form_RPC(67);
			Form_RPC(66);

			FromXY2LonLatTest(470, 3891, 100, lon_lat_mosaic[0], lon_lat_mosaic[1], 67);
			FromlatlonH2xy(lon_lat_mosaic[0], lon_lat_mosaic[1], 100, xy_mosaic[0], xy_mosaic[1], 67);

			cout << lon_lat_mosaic[0] << ",,," << lon_lat_mosaic[1] << endl;

			cout << "From left :(470,3891) -> (" << xy_mosaic[0] << "," << xy_mosaic[1] << ")" << endl;
			int cin_mosaic = 0;
			cin >> cin_mosaic;

			MatchToolType LEFT_INFI_MOSAIC;
			LEFT_INFI_MOSAIC.filepath = LEFT_INFI_MOSAIC.outputpath = "D:\\OpenCVWithCuda\\RPCModel\\ConsoleApplication1\\ConsoleApplication1\\Whole.tif25_.tif";
			MatchToolType RIGHT_INFI_MOSAIC;
			RIGHT_INFI_MOSAIC.filepath = RIGHT_INFI_MOSAIC.outputpath = "D:\\OpenCVWithCuda\\RPCModel\\ConsoleApplication1\\ConsoleApplication1\\Whole.tif26_.tif";

			imgapi1.Mosaic_RigurouModel_Edge(LEFT_INFI_MOSAIC, RIGHT_INFI_MOSAIC, this, 66);
		}



		currentImageID = 62;
        imgapi1.Begin_Match_WorkFlow("C:\\Temp\\match\\left_.tif", "c:\\temp\\match\\right_", this);
        vector<Point2d> _heicalib;
        imgapi1.GetClibContent(_heicalib);
        vector<Point2d> _heicalib_x;
        imgapi1.GetClibContent(_heicalib_x, true);

        calibration.pownumber = 3;
        calibration.CalibrateTime(_heicalib);
        calibration.pownumber = 1;
        calibration.CalibrateTime(_heicalib_x, true);
        vector<TimeCalibStruct> _hei1;
        calibration.IAngleCalibration(_hei1);
        vector<TimeCalibStruct> _hei1_x;
        p_TimeCalib_X = nullptr;
        calibration.IAngleCalibration(_hei1_x, true);
        cout << "畸变数据装载完毕" << endl;
        p_TimeCalib = new vector<TimeCalibStruct>();
        for (int pi = 0; pi < _hei1.size(); pi++)
        {
            p_TimeCalib->push_back(_hei1[pi]);
        }

        p_TimeCalib_X = new vector<TimeCalibStruct>();
        for (int pi = 0; pi < _hei1_x.size(); pi++)
        {
            p_TimeCalib_X->push_back(_hei1_x[pi]);
        }

        cout << "开始校正指向角" << endl;
        base_point = SetPoint(imgframeheight, imgframewidth);
        readcontrolpoint.open("C:\\Temp\\CONTROLPOINTS.json", ios::in);
        imgapi1.Begin_Match_WorkFlow("C:\\Temp\\match\\left_.tif", "c:\\temp\\match\\right_.tiff", this);
        if (!readcontrolpoint.is_open())
        {
            cout << "控制点文件不存在,控制点外畸变参数校正将不被执行,详情请见日志" << endl;
            IAddLog("控制点文件%s不存在，请检查配置，相应的控制点外畸变校正工作将不被畸形，畸变参数将不会与之前发生变化", "C:\\Temp\\CONTROLPOINTS.json");
        }
        do
        {
            char readcontrolLine[1024];
            readcontrolpoint.getline(readcontrolLine, 1024);
            //需要行时 
            ImgPoint newpoints(static_cast<string>(readcontrolLine), asJSon);
            int _lineindex = currentImageID * _singleframeHeight + _singleframeHeight - newpoints.ImgY();
            double time = _simpleTimes[_lineindex].utctime;  //时间的设置 
            newpoints.utcapi = UTCAPI(time, Enum_TimeStandard_ZY03);
            newpoints = newpoints + base_point;
            double xvalue, yvalue, zvalue;
            double heilon, heilat;
            FromXY2LonLatTest(newpoints.ImgX(), newpoints.ImgY(), 100, heilon, heilat, 62);
            pEarth->Geo2Rect(-heilon* CV_PI / 180.0, heilat * CV_PI / 180.0, 100, xvalue, yvalue, zvalue);
            double xposvalue, yposvalue, zposvalue;
            pEarth->Geo2Rect(-newpoints.Lon() * CV_PI / 180.0, newpoints.Lat() * CV_PI / 180.0, 100, xposvalue, yposvalue, zposvalue);

            double test1, test2, test3;
            newpoints.OriX(xvalue);
            newpoints.OriY(yvalue);
            newpoints.OriZ(zvalue);
            newpoints.X(xposvalue); newpoints.Y(yposvalue); newpoints.Z(zposvalue);
            newpoints.Serialize("d:\\testPos\\0410\\");
            double heisqrt = 0;
            heisqrt += pow(newpoints.X() - newpoints.OriX(), 2.0);
            heisqrt += pow(newpoints.Y() - newpoints.OriY(), 2.0);
            heisqrt += pow(newpoints.Z() - newpoints.OriZ(), 2.0);
            heisqrt = sqrt(heisqrt);
            cout << heisqrt << "\t" << heilon << "\t" << heilat << "\t" << newpoints.Lon() << "\t" << newpoints.Lat() << endl;
            if (fabs(newpoints.Lat() - 0) < 0.001)
                continue;
            //controlist.push_back(newpoints);
        } while (readcontrolpoint.peek() != EOF);
        /*
        imgapi1.Begin_Match_WorkFlow("D:\\OpenCVWithCuda\\RPCModel\\ConsoleApplication1\\ConsoleApplication1\\Whole.tif21_.tif",
            "D:\\IRS\\JB11-1_IRS_000116063_003_002_L0\\JB11-1_IRS_000116063_003_004_002_L3C_TOABT.tiff",
            this);
            */

        return;
        vector<string> imagenames;
        ///分离出影像数据
        for (int i = 0; i < framecount; i++)
        {
            cout << "Progress:" << i + 1 << "/" << framecount << "\r";
            char reader[1024];
            sprintf_s(reader, "%dsubframe.tiff", i);
#ifndef _ImageOutPutCheckingStart_Mode_Test
            if (i > framefindstart)
            {
                //unsigned int* imgpointerorder = nullptr;
                //envAPI.IFormImage(imgpointerorder);
                //continue;
                char reader1[1024];
                sprintf_s(reader1, "%dsubframe.tiff", i);
                imagenames.push_back(static_cast<string>(reader1));
            }
#endif
            ImgFacotryAPI imgapi(imgframewidth, imgframeheight, bytecount, reader);
            unsigned int* imgpointer = imgapi.IRow(0);
            if (!envAPI.IFormImage(imgpointer))
            {
                return;
            }
            if (i > framefindstart)
                imgapi.IOutPut(true);
            else
                imgapi.IOutPut(false);
            if (i == framecount - 1)
            {
                ImgFacotryAPI imgapitemp(imgframewidth, imgframeheight, bytecount, reader);
                imgapitemp.Mosaic(imagenames, imgframewidth * imagenames.size(), imgframeheight, false);
                imgapitemp.DetectorHistoMatchCalibration<unsigned short>("c", 480, 10876, 12, GDT_UInt16);
            }
        }
    }
}


void VirtualMachine::FillPointWithModelValue(vector<ImgPoint> pointarray)
{
    int pointsize = pointarray.size();
    for (int i = 0; i < pointsize; i++)
    {

    }
}


ImageModel VirtualMachine::FormImageModel()
{
    if (_myFactory.size() <= 0)
    {
        cout << "辅助数据获取失败,无法创建成像结构，请检查" << endl;
        IAddLog("辅助数据获取失败,无法创建成像结构，请检查");
    }

    int sizeofAux = _myFactory.size();
    /*
    .用Placement new 解决buffer的问题
    问 题描述：用new分配的数组缓冲时，由于调用了默认构造函数，因此执行效率上不佳。
    若没有默认构造函数则会发生编译时错误。如果你想在预分配的内存上创建 对象，
    用缺省的new操作符是行不通的,比如OrbitAPI，
    AttItudeAPI，LineTimeAPI要解决这个问题，可以用placement new构造。它允许你构造一个新对象到预分配的内存上。*/
    //在栈上提前分配(预留)缓存
    char* buffer_OrbitAPI = new char[sizeof(OrbitAPI)* (sizeofAux + 1)];
    char* buffer_AttitudeAPI = new char[sizeof(AttitudeAPI)* (sizeofAux + 1)];
    char* buffer_LineTime = new char[sizeof(LineTimeAPI)* (sizeofAux + 1)];

    OrbitAPI* pOrbit = (OrbitAPI*)buffer_OrbitAPI;
    AttitudeAPI* pAttitude = (AttitudeAPI*)buffer_AttitudeAPI;
    LineTimeAPI* pLineTime = (LineTimeAPI*)buffer_LineTime;

    for (int i = 0; i < sizeofAux; i++)
    {
        ///不理会返回的指针，因为不需要使用它们
        new (pOrbit + i)OrbitAPI(_myFactory[i].Orbit());
        new (pAttitude + i)AttitudeAPI(_myFactory[i].Attitude());
        new (pLineTime + i)LineTimeAPI(_myFactory[i].LineTime());
    }

    new (pOrbit + sizeofAux)OrbitAPI("nullnode");
    new (pAttitude + sizeofAux)AttitudeAPI("nullnode");
    new (pLineTime + sizeofAux)LineTimeAPI("nullnode");

    //if (modelinput != nullptr) modelinput = nullptr;
    ////在此函数不用对这几个缓存进行析构，交给模型类去完成辅助数据的回收了
    char* buffer_ImageModel = new char[sizeof(ImageModel)];
    ImageModel* modelinput = (ImageModel*)buffer_ImageModel;
    new (modelinput)ImageModel(pOrbit, pAttitude, pLineTime);

    ImageModel returnmodel(*modelinput);
    return returnmodel;
}

double VirtualMachine::GetRangeForOneStep(int imageid)
{
    imageid -= 41;
    //currentImageID = imageid;
    if (imageid == -1) imageid = currentImageID;

    double imgx = 0;
    double imgy = _singleframeHeight;
    ///根据行号和景ID确定时间 
    int _lineindex = imageid * _singleframeHeight + _singleframeHeight - imgy;
    int _oriline = imageid * _singleframeHeight;
    double time = _simpleTimes[_lineindex].utctime;  //时间的设置 
    //cout << "whole time is" << _simpleTimes[(imageid + 1) * _singleframeHeight].utctime - _simpleTimes[imageid * _singleframeHeight].utctime << endl;
    //cout << "starttimesi:" << timestartapi.year << "-" << timestartapi.month << "-" << timestartapi.day << "-" << timestartapi.hour << "-" << timestartapi.minute << "-" << timestartapi.second << endl;

    //cout << "Average time is" << _simpleTimes[(imageid + 1) * _singleframeHeight].avgtime<<endl;
    double body2WGS84Rotation[9] = { 0 };
    UTCAPI timecheck(time, Enum_TimeStandard_ZY03);
    double position[3] = { 0 };
    m_pImgModel->IBody2WGS84Test(time, body2WGS84Rotation, position);
    //cout << "body 2 wgs84 finished" << endl;

    double lightvector[3] = { 0 };
    base_point.ILightVector(imgx, imgy, lightvector);
    double earth_LonLatposition[3] = { 0 };
    double earth_RectPosition[3] = { 0 };
    double lookvector[3] = { 0 };

    CSatOrbit helper;
    //helper.invers_matrix(body2WGS84Rotation, 3);
    double rotmirror[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    double rotcamera[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    //{ 0.0141, 89.9914, 90.0111,
    //90.0887, 0.0292, 89.9721,
    //89.9888, 90.0227, 0.0253 };

    double rotaftercamera[9]; m_base.Multi(rotcamera, body2WGS84Rotation, rotaftercamera, 3, 3, 3);
    double rotaftermirror[9]; m_base.Multi(rotmirror, rotaftercamera, rotaftermirror, 3, 3, 3);


    //memcpy(body2WGS84Rotation, rotaftermirror, sizeof(double)* 9);




    double roll = 0;
    double pitch = 0;
    double yaw = 0;// -1.4 + 2.38594403;
    double rotvecvalue[3] = { roll / 180.0 * CV_PI, pitch / 180.0 * CV_PI, yaw / 180.0 * CV_PI };
    Mat rotatevector(3, 1, CV_64F, &rotvecvalue);



    double rotvalue[9];
    CvMat value1 = rotatevector;
    Mat rotationMatrix(3, 3, CV_64F, &rotvalue);
    CvMat value2 = rotationMatrix;
    cvRodrigues2(&value1, &value2);

    double swingrot[9];
    double startpos = (-2.0875 * CV_PI / 180.0);
    double swingoffset = 0;// (0.5* 0.0381971863420548805845321032094 / (9226.0 - / 9226.0 / 9226.0) * CV_PI / 180.0 * imgy*imgy*imgy;
    double wholevalue = +fabs(startpos) * 2 / (1.0006905003015936 * 1) * 0.999811018;
    double imgyshift = 0;
    if (p_TimeCalib != nullptr)
        imgyshift = GetAlongAngle(imgy);

    double imgx_yshift = 0;
    if (p_TimeCalib_X != nullptr)
    {
        imgx_yshift = GetAlongAngle(imgx, true);
    }
    double swingangle1 = startpos + imgy / 10785.0 * wholevalue;
    swingangle1 -= swingoffset;

    double swingan[3] = { -swingangle1 * 2 + (imgyshift - imgx_yshift) / 10785 * wholevalue, 0, 0 };
    double rangelon[2], rangelat[2];
    for (int i = 0; i < 2; i++)
    {
        if (i == 1)
        {
            swingan[0] = -swingan[0];
        }
        //cout << "swing angle is " << swingangle1 << "iMAGEY IS" << imgy * 180.0 / CV_PI << endl;
        //m_base.rot(0, swingan,0, swingrot);
        Mat swingvector(3, 1, CV_64F, &swingan);
        CvMat values1 = swingvector;
        Mat rotationMatrixS(3, 3, CV_64F, &swingrot);
        CvMat values2 = rotationMatrixS;
        cvRodrigues2(&values1, &values2);
        double rotafterswing[9];
        //helper.rot(0, swingangle1, 0, swingrot);
        //m_base.invers_matrix(swingrot, 3);

        double rotafterCompensate[9];
        helper.mult(rotaftermirror, rotvalue, rotafterCompensate, 3, 3, 3);
        helper.mult(rotafterCompensate, swingrot, body2WGS84Rotation, 3, 3, 3);

        //double lightvectorbeforehand[3] = { 0 };
        helper.mult(body2WGS84Rotation, lightvector, lookvector, 3, 3, 1);
        //helper.mult(lightvectorbeforehand, swingrot, lookvector, 1, 3, 3);
        pEarth->IPosition(position, lookvector, 100, earth_RectPosition, earth_LonLatposition);

        rangelon[i] = earth_LonLatposition[1];
        rangelat[i] = earth_LonLatposition[0];
    }
    
    double xrange[2], yrange[2], zrange[2];
    pEarth->Geo2Rect(-rangelon[0], rangelat[0], 100, xrange[0], yrange[0], zrange[0]);
    pEarth->Geo2Rect(-rangelon[1], rangelat[1], 100, xrange[1], yrange[1], zrange[1]);

    return sqrt(pow(xrange[1] - xrange[0], 2.0) + pow(yrange[1] - yrange[0], 2.0) + pow(zrange[1] - zrange[0], 2.0));

    
}

void VirtualMachine::FromXY2LonLat(const double& imgx, const double& imgy, const double& height, double& lon, double& lat, int imageid)
{
    if (imageid == -1) imageid = currentImageID;
    else if (imageid < 0 || imageid >  _myFactory[0].framecount)
    {
        cout << "影像索引号无效，请检查，此次地面坐标求解过程失败" << endl;
        IAddLog("影像索引号无效，请检查，此次地面坐标求解过程失败");
        return;
    }
    ImgPoint imgpoint = base_point;

    ///根据行号和景ID确定时间 
    int _lineindex = imageid * _singleframeHeight + imgy;
    double time = _simpleTimes[_lineindex].utctime;  //时间的设置 
    double body2WGS84Rotation[9] = { 0 };
    double position[3] = { 0 };
    m_pImgModel->IBody2WGS84(time, body2WGS84Rotation, position);
    double lightvector[3] = { 0 };
    imgpoint.ILightVector(imgx, imgy, lightvector);
    double earth_LonLatposition[3] = { 0 };
    double earth_RectPosition[3] = { 0 };
    double lookvector[3] = { 0 };

    CSatOrbit helper;
    helper.invers_matrix(body2WGS84Rotation, 3);
    helper.mult(lightvector, body2WGS84Rotation, lookvector, 1, 3, 3);
    pEarth->IPosition(position, lookvector, height, earth_RectPosition, earth_LonLatposition);
    earth_LonLatposition[0] = earth_LonLatposition[0] / 3.1415926535897932384626433832795 * 180.0;
    earth_LonLatposition[1] = earth_LonLatposition[1] / 3.1415926535897932384626433832795 * 180.0;
}

double VirtualMachine::distanceOnePixel(int indexid, const double& x, const double& y, const double& H, const double& dx, const double& dy)
{
    double lat1, lat2, lon1, lon2;
    FromXY2LonLatTest( x, y, H, lat1, lon1,indexid);
    FromXY2LonLatTest( x + dx, y + dy, H, lat2, lon2,indexid);
    double dis = sqrt((lat1 - lat2)*(lat1 - lat2) + (lon1 - lon2)*(lon1 - lon2)) / sqrt(dx*dx + dy*dy);
    return dis;
}

bool VirtualMachine::FromlatlonH2xy(double lat, double lon, double H, double &x, double &y, int ccdID)
{
	if (rpc_container[ccdID - 41].Is_Inited())
	{
		rpc_container[ccdID - 41].LightIn(lat, lon, H, x, y);

		return true;
	}

	cout << "Lonlat!!!!!" << endl;
    double latp, lonp;
    PredictbasedAffine(ccdID, lat, lon, H, x, y);
    if (x < -20 || x > _singleframeWidth - 1 + 20) return false;
    if (y < -50 || y > _singleframeHeight - 1 + 50) return false;
    //cout << "预测点:" << "(" << x << "," << y << ")" << endl;
    FromXY2LonLatTest( x, y, H, latp, lonp,ccdID);
    double pixelsize = distanceOnePixel(ccdID, x, y, H, 10, 10);
    //cout << "debug:分辨率为:" << pixelsize << endl;
    double e, e0;
    e = sqrt((latp - lat)*(latp - lat) + (lonp - lon)*(lonp - lon)) / pixelsize;


    int numR = 0;
    do {
        numR++;
        e0 = e;
        FromXY2LonLatTest(x, y, H, latp, lonp, ccdID);
		//cout << "x,y" << x << "," << y << endl;
        //cout << "lon,lat 误差" << fabs((lonp - lon) / pixelsize) << "," << fabs((latp - lat) / pixelsize) << endl;
        pixelsize = distanceOnePixel(ccdID, x, y, H, fabs((lonp - lon) / pixelsize), fabs((latp - lat) / pixelsize));
		//cout << "pixelsize is :" << pixelsize << endl;
		double max_step = max(fabs((lonp - lon) / pixelsize), fabs((latp - lat) / pixelsize));
		max_step = max(0.01, max_step);
		if (max_step == 0) break;
		PreciseBasedAffine(ccdID, x, y, lat, lon, H, max_step, max_step);
		//cout << "修正后x,y" << x << "," << y << endl;
		FromXY2LonLatTest(x, y, H, latp, lonp, ccdID);
		//cout << "修正前pixelsize is :" << pixelsize << endl;
        pixelsize = distanceOnePixel(ccdID, x, y, H, fabs((lonp - lon) / pixelsize), fabs((latp - lat) / pixelsize));
		//cout << "修正后pixelsize is :" << pixelsize << endl;
		e = sqrt((latp - lat)*(latp - lat) + (lonp - lon)*(lonp - lon)) / pixelsize;

        if (e<0.001) break;

        if (numR>3)
        {
            return true;
        }

    } while (true);

    return true;
}

bool VirtualMachine::PreciseBasedAffine(int indexid,double &x, double &y, double latitude, double longitude, double H, double dx, double dy)
{

    double s[5], l[5], lat[5], lon[5];

    s[0] = x; l[0] = y;
    s[1] = x + dx; l[1] = y + dy;
    s[2] = x + dx; l[2] = y - dy;
    s[3] = x - dx; l[3] = y + dy;
    s[4] = x - dx; l[4] = y - dy;

    for (int i = 0; i<5; i++)
        FromXY2LonLatTest(s[i], l[i], H, lat[i], lon[i],indexid);

    double LAT_OFF, LONG_OFF, LINE_OFF, SAMP_OFF, LAT_SCALE, LONG_SCALE, LINE_SCALE, SAMP_SCALE;

    orbithelper.Compute_avAnddx(lat, 5, LAT_OFF, LAT_SCALE);
    orbithelper.Compute_avAnddx(lon, 5, LONG_OFF, LONG_SCALE);
    orbithelper.Compute_avAnddx(l, 5, LINE_OFF, LINE_SCALE);
    orbithelper.Compute_avAnddx(s, 5, SAMP_OFF, SAMP_SCALE);

    //¹éÒ»»¯
    for (int i = 0; i<5; i++)
    {
        l[i] = (l[i] - LINE_OFF) / LINE_SCALE;
        s[i] = (s[i] - SAMP_OFF) / SAMP_SCALE;

        lat[i] = (lat[i] - LAT_OFF) / LAT_SCALE;
        lon[i] = (lon[i] - LONG_OFF) / LONG_SCALE;

    }
    //
    //s=ps[0]+ps[1]*lat+ps[2]*lon;
    //l=pl[0]+pl[1]*lat+pl[2]*lon;
    double aas[9], aal[9], a[3], ps[3], pl[3];
    memset(aas, 0, sizeof(double)* 9);
    memset(aal, 0, sizeof(double)* 9);
    memset(ps, 0, sizeof(double)* 3);
    memset(pl, 0, sizeof(double)* 3);
    for (int i = 0; i<5; i++)
    {
        a[0] = 1, a[1] = lat[i], a[2] = lon[i];
        orbithelper.pNormal(a, 3, s[i], aas, ps, 1);
        orbithelper.pNormal(a, 3, l[i], aal, pl, 1);
    }
    orbithelper.Gauss(aas, ps, 3);
    orbithelper.Gauss(aal, pl, 3);

    latitude = (latitude - LAT_OFF) / LAT_SCALE;
    longitude = (longitude - LONG_OFF) / LONG_SCALE;
    x = ps[0] + ps[1] * latitude + ps[2] * longitude;
    y = pl[0] + pl[1] * latitude + pl[2] * longitude;

    y = y*LINE_SCALE + LINE_OFF;
    x = x*SAMP_SCALE + SAMP_OFF;
    return true;
    /*
    double s[5], l[5], lat[5], lon[5];
    s[0] = (_singleframeWidth - 1) / 2; l[0] = (_singleframeHeight - 1) / 2;
    s[1] = 0; l[1] = 0;
    s[2] = _singleframeWidth - 1; l[2] = 0;
    s[3] = _singleframeWidth - 1; l[3] = _singleframeHeight - 1;
    s[4] = 0; l[4] = _singleframeHeight - 1;

    for (int i = 0; i<5; i++)
        FromXY2LonLatTest(s[i], l[i], H, lat[i], lon[i], indexid);

    double LAT_OFF, LONG_OFF, LINE_OFF, SAMP_OFF, LAT_SCALE, LONG_SCALE, LINE_SCALE, SAMP_SCALE;


    orbithelper.Compute_avAnddx(lat, 5, LAT_OFF, LAT_SCALE);
    orbithelper.Compute_avAnddx(lon, 5, LONG_OFF, LONG_SCALE);
    orbithelper.Compute_avAnddx(l, 5, LINE_OFF, LINE_SCALE);
    orbithelper.Compute_avAnddx(s, 5, SAMP_OFF, SAMP_SCALE);

    //¹éÒ»»¯
    for (int i = 0; i<5; i++)
    {
        l[i] = (l[i] - LINE_OFF) / LINE_SCALE;
        s[i] = (s[i] - SAMP_OFF) / SAMP_SCALE;

        lat[i] = (lat[i] - LAT_OFF) / LAT_SCALE;
        lon[i] = (lon[i] - LONG_OFF) / LONG_SCALE;

    }
    //
    //s=ps[0]+ps[1]*lat+ps[2]*lon;
    //l=pl[0]+pl[1]*lat+pl[2]*lon;
    double aas[9], aal[9], a[3], ps[3], pl[3];
    memset(aas, 0, sizeof(double)* 9);
    memset(aal, 0, sizeof(double)* 9);
    memset(ps, 0, sizeof(double)* 3);
    memset(pl, 0, sizeof(double)* 3);
    for (int i = 0; i<5; i++)
    {
        a[0] = 1, a[1] = lat[i], a[2] = lon[i];
        orbithelper.pNormal(a, 3, s[i], aas, ps, 1);
        orbithelper.pNormal(a, 3, l[i], aal, pl, 1);
    }

    orbithelper.Gauss(aas, ps, 3);
    orbithelper.Gauss(aal, pl, 3);

    latitude = (latitude - LAT_OFF) / LAT_SCALE;
    longitude = (longitude - LONG_OFF) / LONG_SCALE;
    x = ps[0] + ps[1] * latitude + ps[2] * longitude;
    y = pl[0] + pl[1] * latitude + pl[2] * longitude;

    y = y*LINE_SCALE + LINE_OFF;
    x = x*SAMP_SCALE + SAMP_OFF;
    return TRUE;
    */
}

bool VirtualMachine::PredictbasedAffine(int indexid,double latitude, double longitude, double H, double &x, double &y)
{

    int line, sample;
    line = _singleframeHeight - 1;
    sample = _singleframeWidth - 1;

    double s[5], l[5], lat[5], lon[5];
    s[0] = sample / 2; l[0] = line / 2;
    s[1] = 0; l[1] = 0;
    s[2] = sample; l[2] = 0;
    s[3] = sample; l[3] = line;
    s[4] = 0; l[4] = line;

    for (int i = 0; i<5; i++)
        FromXY2LonLatTest( s[i], l[i], H, lat[i], lon[i],indexid);

    double LAT_OFF, LONG_OFF, LINE_OFF, SAMP_OFF, LAT_SCALE, LONG_SCALE, LINE_SCALE, SAMP_SCALE;

    orbithelper.Compute_avAnddx(lat, 5, LAT_OFF, LAT_SCALE);
    orbithelper.Compute_avAnddx(lon, 5, LONG_OFF, LONG_SCALE);
    orbithelper.Compute_avAnddx(l, 5, LINE_OFF, LINE_SCALE);
    orbithelper.Compute_avAnddx(s, 5, SAMP_OFF, SAMP_SCALE);

    //¹éÒ»»¯
    for (int i = 0; i<5; i++)
    {
        l[i] = (l[i] - LINE_OFF) / LINE_SCALE;
        s[i] = (s[i] - SAMP_OFF) / SAMP_SCALE;

        lat[i] = (lat[i] - LAT_OFF) / LAT_SCALE;
        lon[i] = (lon[i] - LONG_OFF) / LONG_SCALE;

    }
    //
    //s=ps[0]+ps[1]*lat+ps[2]*lon;
    //l=pl[0]+pl[1]*lat+pl[2]*lon;
    double aas[9], aal[9], a[3], ps[3], pl[3];
    memset(aas, 0, sizeof(double)* 9);
    memset(aal, 0, sizeof(double)* 9);
    memset(ps, 0, sizeof(double)* 3);
    memset(pl, 0, sizeof(double)* 3);
    for (int i = 0; i<5; i++)
    {
        a[0] = 1, a[1] = lat[i], a[2] = lon[i];
        orbithelper.pNormal(a, 3, s[i], aas, ps, 1);
        orbithelper.pNormal(a, 3, l[i], aal, pl, 1);
    }
    orbithelper.Gauss(aas, ps, 3);
    orbithelper.Gauss(aal, pl, 3);

    latitude = (latitude - LAT_OFF) / LAT_SCALE;
    longitude = (longitude - LONG_OFF) / LONG_SCALE;
    x = ps[0] + ps[1] * latitude + ps[2] * longitude;
    y = pl[0] + pl[1] * latitude + pl[2] * longitude;

    y = y*LINE_SCALE + LINE_OFF;
    x = x*SAMP_SCALE + SAMP_OFF;
    return TRUE;
}


bool VirtualMachine::GetNeightbor(int indexid,vector<Point2d>& _LEFTMATCH,vector<Point2d>& _RIGHTMATCH)
{
    double x_index, y_index;
    y_index = 100;
    do
    {
        y_index += 3;
        x_index = 15;
        do
        {
            
            double loninfo,  latinfo;
            FromXY2LonLatTest(x_index, y_index, 100, loninfo, latinfo, indexid + 1);
            double xinfo, yinfo;
            FromlatlonH2xy(latinfo, loninfo, 100, xinfo, yinfo, indexid);
            if (xinfo > _singleframeWidth - 1 || xinfo < 0)
                break;
            if (yinfo > _singleframeHeight - 1 || yinfo < 0)
                break;

            if (x_index > _singleframeWidth - 1)
                break;
            
            _RIGHTMATCH.push_back(Point2d(x_index,y_index));
            _LEFTMATCH.push_back(Point2d(xinfo,yinfo));
            x_index += 3;
        } while (true);
    } while (y_index + 3 < _singleframeHeight -1);

    ImgFacotryAPI imgapitemp(480, 10786, 12, "hei");
    MatchToolType left_tool, right_tool;
    
    int fileindex = indexid - 41;
    char readerl[1024];
    char readerr[1024];
    sprintf_s(readerl, "D:\\OpenCVWithCuda\\RPCModel\\ConsoleApplication1\\ConsoleApplication1\\Whole.tif%d_.tif", fileindex);
    left_tool.filepath = readerl;
    left_tool.outputpath = readerr;

    sprintf_s(readerr, "D:\\OpenCVWithCuda\\RPCModel\\ConsoleApplication1\\ConsoleApplication1\\Whole.tif%d_.tif", fileindex + 1);
    right_tool.filepath = readerr;
    right_tool.outputpath = readerr;

    left_tool._ptool = static_cast<void*>(this);
    right_tool._ptool = static_cast<void*>(this);

    imgapitemp.GetNeightContent(_LEFTMATCH, _RIGHTMATCH);

    imgapitemp.LSM_Match(left_tool, right_tool, 19);
}


void VirtualMachine::FromXY2LonLatTest(const double& imgx, const double& imgy, const double& height, double& lon, double& lat, int imageid)
{
    imageid -= 41;
	//safestring::showMemoryInfo();
    //currentImageID = imageid;
    if (imageid == -1) imageid = currentImageID;

	if (rpc_container[imageid].Is_Inited())
	{
		rpc_container[imageid].LightOut(imgx, imgy, height, lon, lat);
		//cout << "lon:" << lon << "lat: " << lat << endl;
		return;
	}
	
    
    ///根据行号和景ID确定时间 
    int _lineindex = imageid * _singleframeHeight + _singleframeHeight- 1 - imgy;
    int _oriline = imageid * _singleframeHeight;
	double time;

	int firstindex = imageid * _singleframeHeight + _singleframeHeight - 1;
	///第一行第二行时间之差
	double time_step = _simpleTimes[firstindex].utctime - _simpleTimes[firstindex - 1].utctime;
	if (imgy < 0)
	{
		time = _simpleTimes[firstindex].utctime + imgy * (-time_step);
	}
	else if (imgy > _singleframeHeight -1)
	{
		int lastindex = imageid * _singleframeHeight;
		double time_step = _simpleTimes[lastindex].utctime - _simpleTimes[lastindex + 1].utctime;
		time = _simpleTimes[lastindex].utctime + (imgy - _singleframeHeight + 1) * time_step;
	}
	else
	 time = _simpleTimes[_lineindex].utctime;  //时间的设置 
    double addon1 =  10 - fabs(imgy - 5392.5) / 5392.5 * 10.0;
    time += addon1 / 480.0 *   10786 * 3.5122304903644342682694142339248e-5;
	if (imageid == 24)
		time -= 40 / 480.0 * 10786 * 3.5122304903644342682694142339248e-5;
	if (imageid == 26)
		time += 43 / 480.0 * 10786 * 3.5122304903644342682694142339248e-5;

    //time -= 10 * 3.5122304903644342682694142339248e-5;
    //cout.precision(20);
    //cout << "----current utc time is: " << time<<"----";

    //cout << "whole time is" << _simpleTimes[(imageid + 1) * _singleframeHeight].utctime - _simpleTimes[imageid * _singleframeHeight].utctime << endl;
    //cout << "starttimesi:" << timestartapi.year << "-" << timestartapi.month << "-" << timestartapi.day << "-" << timestartapi.hour << "-" << timestartapi.minute << "-" << timestartapi.second << endl;

    //cout << "Average time is" << _simpleTimes[(imageid + 1) * _singleframeHeight].avgtime<<endl;
    double body2WGS84Rotation[9] = { 0 };
   // UTCAPI timecheck(time, Enum_TimeStandard_ZY03);
    double position[3] = { 0 };
	//safestring::showMemoryInfo();

    m_pImgModel->IBody2WGS84Test(time, body2WGS84Rotation, position);
	//safestring::showMemoryInfo();

    //cout << "body 2 wgs84 finished" << endl;

    double lightvector[3] = { 0 };
    base_point.ILightVector(imgx, imgy, lightvector);
    double earth_LonLatposition[3] = { 0 };
    double earth_RectPosition[3] = { 0 };
    double lookvector[3] = { 0 };

    //CSatOrbit helper;
    //helper.invers_matrix(body2WGS84Rotation, 3);
    double rotmirror[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    double rotcamera[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    //{ 0.0141, 89.9914, 90.0111,
    //90.0887, 0.0292, 89.9721,
    //89.9888, 90.0227, 0.0253 };

    double rotaftercamera[9]; m_base.Multi(rotcamera, body2WGS84Rotation, rotaftercamera, 3, 3, 3);
    double rotaftermirror[9]; m_base.Multi(rotmirror, rotaftercamera, rotaftermirror, 3, 3, 3);


    //memcpy(body2WGS84Rotation, rotaftermirror, sizeof(double)* 9);




    double roll = 0;
    double pitch = 0;
    double yaw = 0;// -1.4 + 2.38594403;
    double rotvecvalue[3] = { roll / 180.0 * CV_PI, pitch / 180.0 * CV_PI, yaw / 180.0 * CV_PI };
    Mat rotatevector(3, 1, CV_64F, &rotvecvalue);



    double rotvalue[9];
    CvMat value1 = rotatevector;

    Mat rotationMatrix(3, 3, CV_64F, &rotvalue);
    CvMat value2 = rotationMatrix;
    cvRodrigues2(&value1, &value2);

    double swingrot[9];
    double startpos = (-2.0875 * CV_PI / 180.0);
    double swingoffset = 0;// (0.5* 0.0381971863420548805845321032094 / (9226.0 - / 9226.0 / 9226.0) * CV_PI / 180.0 * imgy*imgy*imgy;
    double wholevalue = +fabs(startpos) * 2 * 1.00549388739566009032939488764 * 1.000250944846080260950021170825 / 1.0057495352300516595331044322155 * 1.000931158 / 1.0001140023647047061221734279613;
        ;
    double imgyshift = 0;


    double imgx_yshift = 0;
	
    if (p_TimeCalib_X != nullptr)
    {
        imgx_yshift = GetAlongAngle(imgx, true);
    }


    double valueget = 0;
	/*
    if (_gofuckyourselfSwingError.size() > 0)
    {
        //cout << "hei" << endl;
        if (imgy > _gofuckyourselfSwingError[0]._start && imgy < endswingvalue)
        {
            

            if (imgy < _gofuckyourselfSwingError[_gofuckyourselfSwingError.size() - 1]._start)
            {
                double heibegin = _gofuckyourselfSwingError[_gofuckyourselfSwingError.size() - 1]._start;
                int step = (heibegin - imgy) / 5.0;
                

                for (int i = 0; i < step; i++)
                {
                    valueget += 5.0 * _gofuckyourselfSwingError[_gofuckyourselfSwingError.size() - 1 - 1 - i]._avg;
                }

                valueget += _gofuckyourselfSwingError[_gofuckyourselfSwingError.size() - 2 - step]._avg * fabs(heibegin - 5.0 * step - imgy);

                valueget += _gofuckyourselfSwingError[_gofuckyourselfSwingError.size() - 1]._avg * (endswingvalue - heibegin);
            }

            else
            {
                valueget += _gofuckyourselfSwingError[_gofuckyourselfSwingError.size() - 1]._avg * fabs(endswingvalue - imgy);
            }
        }
    }
	*/
    double swingangle1 = 0;
	
    if (true/*fabs(valueget - 0) < 0.0001*/)
    {
        double addon = 6 - fabs(imgy - 5392.5) / 5392.5 * 6;
		if (imageid == 26)
			addon -= 5;
        swingangle1 = (startpos + (10785 - imgy ) / 10785.0 * wholevalue) * 2 + (-5.0 - addon)/ 10785.0 * wholevalue * 2;
		
	}

    else
    {
        //cout << "hei" << endl;
        swingangle1 = startpos * 2 + (10785 - endswingvalue) / 10785 * wholevalue * 2 + valueget / 10785 * wholevalue * 2;
    }
        //swingangle1 -= swingoffset;
    //double realswan = -swingangle1 * 2 + 0;
    double swingan[3] = { swingangle1, 0, 0 };
    //cout << "swing angle is " << swingangle1 << "iMAGEY IS" << imgy * 180.0 / CV_PI << endl;
    //m_base.rot(0, swingan,0, swingrot);
    Mat swingvector(3, 1, CV_64F, &swingan);
    CvMat values1 = swingvector;
    Mat rotationMatrixS(3, 3, CV_64F, &swingrot);
    CvMat values2 = rotationMatrixS;
    cvRodrigues2(&values1, &values2);
    double rotafterswing[9];
    //helper.rot(0, swingangle1, 0, swingrot);
    //m_base.invers_matrix(swingrot, 3);

    double rotafterCompensate[9];
	
    orbithelper.mult(rotaftermirror, rotvalue, rotafterCompensate, 3, 3, 3);
	orbithelper.mult(rotafterCompensate, swingrot, body2WGS84Rotation, 3, 3, 3);

    //double lightvectorbeforehand[3] = { 0 };
	orbithelper.mult(body2WGS84Rotation, lightvector, lookvector, 3, 3, 1);
    //helper.mult(lightvectorbeforehand, swingrot, lookvector, 1, 3, 3);
	//safestring::showMemoryInfo();

    pEarth->IPosition(position, lookvector, height, earth_RectPosition, earth_LonLatposition);
	//safestring::showMemoryInfo();

    earth_LonLatposition[0] = earth_LonLatposition[0] / 3.1415926535897932384626433832795 * 180.0;
    earth_LonLatposition[1] = earth_LonLatposition[1] / 3.1415926535897932384626433832795 * 180.0;


    lon = earth_LonLatposition[1];
    lat = earth_LonLatposition[0];

	rotatevector.~Mat();
	swingvector.~Mat();
	rotationMatrixS.~Mat();
}

void VirtualMachine::FromXY2LonLatUnitTest(const double& imgx, const double& imgy, const double& height, double& lon, double& lat, int imageid, string inputer)
{
	int debug_input;
    imageid -= 62;
    ImgPoint imgpoint = base_point;
    double body2WGS84Rotation[9] = { 0 };
    double position[3] = { 0 };

	cin >> debug_input;
    m_pImgModel->IBody2WGS84(inputer, body2WGS84Rotation, position);
	cin >> debug_input;
    double lightvector[3] = { 0 };
    imgpoint.ILightVector(imgx, imgy, lightvector);
    double earth_LonLatposition[3] = { 0 };
    double earth_RectPosition[3] = { 0 };
    double lookvector[3] = { 0 };


    if (!pEarth->IPosition(position, lookvector, height, earth_RectPosition, earth_LonLatposition))
    {
        cout << "定位测试失败，请检查" << endl;
        IAddLog("定位测试失败，请检查");
    }
    earth_LonLatposition[0] = earth_LonLatposition[0] / 3.1415926535897932384626433832795 * 180.0;
    earth_LonLatposition[1] = earth_LonLatposition[1] / 3.1415926535897932384626433832795 * 180.0;
}

void VirtualMachine::SetEarthModel()
{
    if (pEarth != nullptr) pEarth = nullptr;
    pEarth = new EarthModel(6378137.0, 6356752.314, 0.08181919, 0.006694379, 0.0, 0.0, 0.0);
}

void VirtualMachine::SetSimpleTime(vector<int> frameheight)
{
    int auxsize = _myFactory.size();
    if (auxsize <= 0)
    {
        return;
    }

    int previousLine = 10;
    vector<SimpleTime> _tempsimple;
    int ncountsize = 1;
    for (int i = 0; i < auxsize; i++)
    {
        SimpleTime _nowtime;
        _nowtime.utctime = _myFactory[i].LineTime().IGetTime();
        _nowtime.actLine = _myFactory[i].LineTime().ACTLine() * 0.5;
        _nowtime.lineid = _nowtime.actLine + i * _singleframeHeight;
        _nowtime.frameid = i + 41;

        _myFactory[i].LineTime().ACTLine(_nowtime.actLine);
        _tempsimple.push_back(_nowtime);
    }



    for (int i = 0; i < auxsize; i++)
    {
        vector<SimpleTime> _temparray;
        if (i == 0)
            _temparray.push_back(_tempsimple[i + 1]);
        else
            _temparray.push_back(_tempsimple[i - 1]);
        _temparray.push_back(_tempsimple[i]);

        DetailedLineTimeInfo(_temparray, i + 41);
    }

    cout << "行时计算完毕" << endl;
}


void VirtualMachine::DetailedLineTimeInfo(vector<SimpleTime> _info, int frameagg)
{
    ///第一幅影像
    double timesegment = 0;
    //cout << _info[0].frameid << "," << _info[1].frameid << endl;
    if (_info[0].frameid < _info[1].frameid)
    {
        timesegment = (_info[1].utctime - _info[0].utctime) / (10786 - _info[0].actLine + _info[1].actLine);
    }
    else{
        timesegment = (_info[0].utctime - _info[1].utctime) / (10786 - _info[1].actLine + _info[0].actLine);
    }

    timesegment = 3.5122304903644342682694142339248e-5 / 1;

    ///设置摆扫成像时间
    perT = timesegment;
    for (int i = 0; i < _singleframeHeight; i++)
    {
        SimpleTime timenow;
        timenow.actLine = i;
        timenow.avgtime = timesegment;
        timenow.frameid = frameagg;
        timenow.lineid = (frameagg - 41) * _singleframeHeight + i;
        timenow.utctime = _info[1].utctime + (i - _info[1].actLine) * timesegment;
        _simpleTimes.push_back(timenow);
    }
}

void VirtualMachine::Calibration()
{

}

void VirtualMachine::Form_RPC(int index_id)
{
	char reader[1024]; 
	sprintf_s(reader, "c:\\Temp\\match\\rpc_model_%d.txt", index_id);
	ifstream ifs(reader);
	if (ifs.is_open())
	{
		ifs.close();
		cout << index_id - 41<<"rpc文件存在!" << endl;
		rpc_container[index_id - 41].IFile(reader,true,false);
		return;
	}
	ifs.close();

	cout << "rpc文件不存在，重新生成中" << endl;
	vector<FPoint> fpoint_container;
	int ncounter_v = 0;
	for (int t_z = 20; t_z < 120; t_z += 40)
	for (int t_y = 1; t_y < _singleframeHeight - 1; t_y += 2)
	{
		for (int t_x = 1; t_x < _singleframeWidth - 1; t_x += 2)
		{
			ncounter_v++;
		}
	}

	fpoint_container.reserve(ncounter_v);
	fpoint_container.resize(ncounter_v);
	int sizeofall = 0;
	for (int t_z = 20; t_z < 120; t_z += 40)
	for (int t_y = 1; t_y < _singleframeHeight - 1; t_y += 2)
	{
		cout << "t_z:" << t_z << ":" << t_y << "\\" << _singleframeHeight - 1 << "\r";
		for (int t_x = 1; t_x < _singleframeWidth - 1; t_x += 2, sizeofall++)
		{
			FPoint point_now;
			point_now.Pixel.x = t_x;
			point_now.Pixel.y = t_y;
			point_now.ground.h = t_z;
			point_now.ground.X = point_now.ground.Y = point_now.ground.Z = -99;
			FromXY2LonLatTest(point_now.Pixel.x, point_now.Pixel.y, t_z, point_now.ground.lon, point_now.ground.lat, index_id);
			fpoint_container[sizeofall] = point_now;
		}
	}

	RPCEvaluationHelper rpc;
	rpc.GetCommander(fpoint_container);
	rpc.ISWitch();
	rpc.IConductor(static_cast<RPCalculator>(15));
	int* marks = new int[80];
	memset(marks, 0, sizeof(int)* 80);
	for (int j = 0; j < 80; j++)
	{
		marks[j] = 1;
	}
	//rpc.ISelection(marks);
	rpc.ISWitch();

	char file_reader[1024];
	sprintf_s(file_reader, "c:\\Temp\\match\\rpc_calced_%d_rpc.txt", index_id);
	rpc.FileName(file_reader);
	rpc.IReport(true);

	rpc_container[index_id - 41] = rpc;

	cout << "第" << index_id - 41 << "景影像RPC模型生成" << endl;
}