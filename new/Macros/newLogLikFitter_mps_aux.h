#ifndef NEWLOGLIKFITTER_MPS_AUX_H
#define NEWLOGLIKFITTER_MPS_AUX_H


///////////////////////////////////////////////////////////////////////////////
// MIN POINT MARKER CREATION HELPER
///////////////////////////////////////////////////////////////////////////////

void min_point_marker_helper
(
    TMarker *&mark_min_point_sysx_l,
    TMarker *&mark_min_point_sysx_h,
    TLine *&line_min_point_sysx_l,
    TLine *&line_min_point_sysx_h,
    double min_point[2],
    double min_point_fake_data[2],
    double min_point_sysx_l[2],
    double min_point_sysx_h[2],
    const bool ENABLE_MIN_POINT_SYSx,
    Int_t markstyle,
    Color_t markcolor,
    Double_t marksize,
    Color_t linecolor
)
{


    if(ENABLE_MIN_POINT_SYSx == true)
    {
        double x_orig = 0.0;
        double y_orig = 0.0;
        if(g_mode_fake_data == true)
        {
            x_orig = min_point_fake_data[0];
            y_orig = min_point_fake_data[1];
        }
        else if(g_mode_fake_data == false)
        {
            x_orig = min_point[0];
            y_orig = min_point[1];
        }

        if((min_point_sysx_l[0] != 0.0) &&
           (min_point_sysx_l[1] != 0.0))
        {
            double x_min_point = min_point_sysx_l[0] - min_point_fake_data[0] + x_orig;
            double y_min_point = min_point_sysx_l[1] - min_point_fake_data[1] + y_orig;

            //std::cout << "DRAW MARK SYSxL" << std::endl;
            std::cout << "SYSxL: " << min_point_sysx_l[0] << " " << min_point_sysx_l[1] << std::endl;
            std::cout << "MIN_POINT: " << min_point[0] << " " << min_point[1] << std::endl;
            std::cout << "MIN_POINT FAKE DATA: " << min_point_fake_data[0] << " " << min_point_fake_data[1] << std::endl;
            std::cout << "DRAWING: " << x_min_point << " " << y_min_point << std::endl;

            mark_min_point_sysx_l = new TMarker(x_min_point, y_min_point, markstyle);
            mark_min_point_sysx_l->SetMarkerColorAlpha(markcolor, 1.0);
            mark_min_point_sysx_l->SetMarkerSize(marksize);

            line_min_point_sysx_l = new TLine(x_orig, y_orig, x_min_point, y_min_point);
            line_min_point_sysx_l->SetLineColorAlpha(linecolor, 0.5);
            line_min_point_sysx_l->SetLineWidth(1);
        }

        if((min_point_sysx_h[0] != 0.0) &&
           (min_point_sysx_h[1] != 0.0))
        {
            double x_min_point = min_point_sysx_h[0] - min_point_fake_data[0] + x_orig;
            double y_min_point = min_point_sysx_h[1] - min_point_fake_data[1] + y_orig;

            //std::cout << "DRAW MARK SYSxH" << std::endl;
            std::cout << "SYSxH: " << min_point_sysx_h[0] << " " << min_point_sysx_h[1] << std::endl;
            std::cout << "MIN_POINT: " << min_point[0] << " " << min_point[1] << std::endl;
            std::cout << "MIN_POINT FAKE DATA: " << min_point_fake_data[0] << " " << min_point_fake_data[1] << std::endl;
            std::cout << "DRAWING: " << x_min_point << " " << y_min_point << std::endl;
            mark_min_point_sysx_h = new TMarker(x_min_point, y_min_point, markstyle);
            mark_min_point_sysx_h->SetMarkerColorAlpha(markcolor, 1.0);
            mark_min_point_sysx_h->SetMarkerSize(marksize);

            line_min_point_sysx_h = new TLine(x_orig, y_orig, x_min_point, y_min_point);
            line_min_point_sysx_h->SetLineColorAlpha(linecolor, 0.5);
            line_min_point_sysx_h->SetLineWidth(1);
        }
    }

}


#endif // NEWLOGLIKFITTER_MPS_AUX_H
