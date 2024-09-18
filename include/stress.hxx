#pragma once

struct StressTensor
{
    bool sigma_xx_known;
    bool sigma_xy_known;
    bool sigma_xz_known;
    bool sigma_yy_known;
    bool sigma_yz_known;
    bool sigma_zz_known;

    double sigma_xx;
    double sigma_xy;
    double sigma_xz;
    double sigma_yy;
    double sigma_yz;
    double sigma_zz;
};
