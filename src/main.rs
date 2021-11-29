#![allow(non_snake_case)]
use std::{io::Write, num::NonZeroI128, ops::Add, slice::Windows};
use plotters::prelude::*;

#[allow(non_upper_case_globals)]
const delta :f64 = 5.0e-5;
const ERROR :f64 = 1.0e-3;
const N:u64 = 512 * 2;

enum NAGARE {
    MthL1,
    MthG1MeL1,
    suehiro,
    deguti,
    ka,
    saiteki,
    husoku,
    humei,
}

#[derive(Default)]
struct Nozzle{
    Di : f64,
    Ai : f64,

    Dth : f64,
    Ath : f64,
    xth : f64,

    De : f64,
    Ae : f64,
    xe : f64,
}

#[derive(Clone, Copy, Default)]
struct quantity{
    x : f64,
    Mx : f64,
    rx :f64,
}

impl Nozzle{
    fn new(Di : f64, Dth:f64, xth : f64, De:f64, xe:f64) -> Nozzle {
        let mut nozzle = Nozzle {Di, Dth, xth, De, xe, ..Default::default()};
        nozzle.init();
        return nozzle;
    }
    fn new_by_AeAc(AeAc:f64) -> Nozzle {
        let mut nozzle = Nozzle {Di:20e-3, Dth:10e-3, xth: 15.0e-3, De:AeAc.powf(0.5) * 10e-3, xe: 100.0e-3, ..Default::default()};
        nozzle.init();
        return nozzle;
    }

    //先細部
    fn calc_AxC(&self, x: f64) -> f64 {
        return std::f64::consts::PI / 4.0 * (self.calc_DxC(x)).powf(2.0); 
    }
    fn calc_DxC(&self, x:f64) -> f64{
        let Dx = (self.Dth - self.Di) * x / self.xth + self.Di;
        return Dx;
    }
    //末広部
    fn calc_AxD(&self, x: f64) -> f64{
        return std::f64::consts::PI / 4.0 * self.calc_DxD(x).powf(2.0); 
    }
    fn calc_DxD(&self, x: f64) -> f64{
        let x_xth = x - self.xth; 
        let Dx = (self.De - self.Dth) * x_xth / (self.xe - self.xth) + self.Dth;
        return Dx; 
    }

    fn init(&mut self){
        self.Ai = std::f64::consts::PI / 4.0 * self.Di.powf(2.0);
        self.Ae = std::f64::consts::PI / 4.0 * self.De.powf(2.0); 
        self.Ath = std::f64::consts::PI / 4.0 * self.Dth.powf(2.0); 
        println!("Ae = {:.2}[mm2]   Ath = {:.2}[mm2]", self.Ae*1.0e6, self.Ath*1.0e6);
    }
}
fn calc_pP(kappa : f64, M : f64) -> f64 {
    return ( 1.0 + (kappa - 1.0)/2.0 * M.powf(2.0) ).powf( -kappa / (kappa - 1.0) );
}
fn calc_Pp(kappa : f64, M : f64) -> f64 {
    return ( 1.0 + (kappa - 1.0)/2.0 * M.powf(2.0) ).powf( kappa / (kappa - 1.0) );
}
fn calc_M2_from_M1_by_AA(kappa: f64, M1: f64, A2A1: f64, isSuper: bool) -> f64 {
    let mut M2 :f64 = 1.0;
    let km1 = kappa - 1.0;
    let kp1 = kappa + 1.0;
    let pwr = kp1 / ( 2.0 * km1 );
    if M1 == 1.0{
        let kmkp = km1 / kp1;
        let kkp1 = 2.0 / kp1;
        while (
                A2A1-
                //M1 / M2 * ( ( (kappa - 1.0) * M2.powf(2.0) + 2.0 )/( (kappa - 1.0) * M1.powf(2.0) + 2.0 ) ).powf( (kappa + 1.0)/(2.0 * (kappa - 1.0))
                1.0 / M2 * ( kmkp * M2.powf(2.0) + kkp1 ).powf(pwr)
            ) > ERROR
        {
            if isSuper{
                M2 += delta;
            } else {
                M2 -= delta;
            }
        }
    }
    else {
        while (
                A2A1 -
                //M1 / M2 * ( ( (kappa - 1.0) * M2.powf(2.0) + 2.0 )/( (kappa - 1.0) * M1.powf(2.0) + 2.0 ) ).powf( (kappa + 1.0)/(2.0 * (kappa - 1.0))
                M1 / M2 * ( ( km1 * M2.powf(2.0) + 2.0 )/( km1 * M1.powf(2.0) + 2.0 ) ).powf(pwr)
            ) > ERROR
        {
            if isSuper{
                M2 += delta;
            } else {
                M2 -= delta;
            }
        }
    }
    //println!("M2 = {:.3}", M2);
    return M2;
}

fn calc_M_before_shock(kappa: f64, PeP0: f64) -> f64 {
    //(8.24)
    let mut M1 : f64 = 1.0;
    while (
            ((kappa + 1.0) * M1.powf(2.0) / ((kappa - 1.0) * M1.powf(2.0) + 2.0)).powf(kappa / (kappa - 1.0))
          * ((kappa + 1.0)/(2.0 * kappa *  M1.powf(2.0) - (kappa - 1.0) )).powf( 1.0/(kappa -1.0) )
          - PeP0 ).abs() > ERROR
    {
        M1 += 0.01;   
    }
    println!("上流M1 = {:.3}", M1);
    return M1;
}
fn calc_M_after_shock(kappa: f64, M1: f64) -> f64 {
    //(8.18)
    let M2 = ( ((kappa -1.0) * M1.powf(2.0) + 2.0) / (2.0 * kappa * M1.powf(2.0) - (kappa - 1.0)) ).powf(0.5);
    println!("下流M2 = {:.3}", M2);
    return M2;   
}
fn calc_AAc(kappa : f64, M: f64) -> f64 {
    let AAc = 1.0 / M * ( ((kappa - 1.0)*M.powf(2.0) + 2.0)/ (kappa + 1.0)).powf( (kappa + 1.0)/(2.0 * (kappa - 1.0)) );
    println!("A/Ac = {:.3}", AAc);
    return AAc;
}
fn siki_8_25(kappa : f64, M1: f64) -> f64 {
    //p2: 衝撃は直後の静圧 p0:ノズル上流全圧
    let p2p0 = (2.0 * kappa * M1.powf(2.0) - (kappa - 1.0))/(kappa + 1.0)*(1.0 + (kappa - 1.0)/2.0 * M1.powf(2.0)).powf(-(kappa/ (kappa - 1.0)));
    return p2p0;
}

fn calc_pcp0(kappa : f64) -> f64 {
    let pcp0 = (2.0 / (kappa + 1.0)).powf(kappa / (kappa - 1.0));
    println!("p*/p0 = {:.3}", pcp0);
    return pcp0;
}

fn check_nagare_keitai(nozzle : &Nozzle, kappa : f64, PbP0: f64) -> NAGARE {
    let AeAth = &nozzle.Ae / &nozzle.Ath;
    println!("Ae/Ath = {:.3}", AeAth);

    let pcP0 = calc_pcp0(kappa);

    let Mde = calc_M2_from_M1_by_AA(kappa, 1.0, AeAth, false);
    let pdP0 = calc_pP(kappa, Mde); //チョーク　下流が亜音速
   
    let Mje = calc_M2_from_M1_by_AA(kappa, 1.0, AeAth, true);
    let pjP0 = calc_pP(kappa, Mje); //チョーク　適正膨張
    
    let phP0 = siki_8_25(kappa, Mje);//チョーク　出口で垂直衝撃波

    println!("pd/P0 = {:.3}  pj/P0 = {:.3}    ph/P0 = {:.3}",pdP0, pjP0, phP0);

    if pdP0 <= PbP0 && PbP0 < 1.0 {
        //亜音速流
        if pdP0 < PbP0 && PbP0 < 1.0 {
            println!("全域で亜音速");
            return NAGARE::MthL1;
        }
        else {
            println!("スロートで音速、その他で亜音速");
            return NAGARE::MthG1MeL1;
        }
    }
    else{ 
        //スロートで音速
        if phP0 < PbP0 && PbP0 < pdP0 {
            println!("非等エントロピー");
            println!("ノズルの末広部に垂直衝撃波");

            return NAGARE::suehiro;
        }
        else if phP0 == PbP0 { 
            println!("出口に垂直衝撃波");

            return NAGARE::deguti;
        }
        else if pjP0 < PbP0 && PbP0 < phP0 {
            println!("過膨張噴流");
            return NAGARE::ka;
        }
        else if pjP0 == PbP0 {
            println!("適正膨張");
            return NAGARE::saiteki;
        } 
        else if PbP0 < pjP0 {
            println!("不足膨張");
            return NAGARE::husoku;       
        }
    }
    return NAGARE::humei;
}

fn calc_Mx(nozzle : &Nozzle, kappa : f64, PbP0: f64, nagare: NAGARE) -> Vec<quantity>{
    let mut Quans : Vec<quantity> = Vec::new();
    let dx = &nozzle.xe / N as f64;
    let mut M2 = 0.0;
    let mut AmAc = 0.0;

    if let NAGARE::suehiro = nagare {
        //衝撃波がある場合
        let M1 = calc_M_before_shock(kappa, PbP0);
        M2 = calc_M_after_shock(kappa, M1);
        AmAc = calc_AAc(kappa, M1);//衝撃波発生位置の断面積比
    }

    for i in 1..N {
        let x = dx * i as f64;
        let mut Mx = 0.0;
        let mut rx = 0.0;
        if x <= nozzle.xth {
            //先細部の計算
            // 状態2:x     状態1:th
            let AxAth = nozzle.calc_AxC(x) / nozzle.Ath;
            rx = nozzle.calc_DxC(x) / 2.0;
            Mx = calc_M2_from_M1_by_AA(kappa, 1.0, AxAth, false);
        }
        else {
            let AxAth = nozzle.calc_AxD(x) / nozzle.Ath;
            rx = nozzle.calc_DxD(x) / 2.0;
            match nagare {
                NAGARE::suehiro => {
                    //末広部の計算
                    //非等エントロピー
                    if AxAth < AmAc {
                        //衝撃波の上流 
                        Mx = calc_M2_from_M1_by_AA(kappa, 1.0, AxAth, true);
                    }
                    else {
                        //衝撃波の下流
                        let AxAm = AxAth / AmAc;
                        Mx = calc_M2_from_M1_by_AA(kappa, M2, AxAm, false);            
                    }
                },
                NAGARE::MthG1MeL1 | NAGARE::MthL1 => {
                    //等エントロピー、亜音速
                    Mx = calc_M2_from_M1_by_AA(kappa, 1.0, AxAth, false);
                },
                _ => {
                    //等エントロピー、超音速
                    Mx = calc_M2_from_M1_by_AA(kappa, 1.0, AxAth, true);
                }
            }
        }
        
        if i == N-2{
            println!("Me = {:.3}", Mx);
        }
        let res = quantity {x, Mx, rx};
        Quans.push(res);
    }
    return Quans;
}
fn main() -> Result<(), Box<std::error::Error>>{
    //let mut nozzle = Nozzle::new(20.0e-3, 10.0e-3, 10.0e-3, 14.1421356e-3, 40.0e-3);
    let mut nozzle = Nozzle::new_by_AeAc(4.0);

    let kappa : f64 = 1.40;
    let Mw = 28.8;
    let RR = 8.3143e3;//J/kmolK
    let R = RR/Mw;
    println!("R = {:.3}[kJ/kgK]", R/1000.0);

    //let mut P0 = 102.0e3;//全域亜音速
    //let P0 = 150.0e3;//衝撃波
    let P0 = 1000e3;
    let Pa = 101.325e3;//Pa
    let mut PbP0 = Pa / P0;
    println!("PbP0 = {:.3}", PbP0);

    let nagare = check_nagare_keitai(&nozzle, kappa, PbP0);
    let mut quans = calc_Mx(&nozzle, kappa, PbP0, nagare);

    /*
    let mut i = 0;
    while PbP0 > 0.4  {
        PbP0 = Pa / P0;

        let nagare = check_nagare_keitai(&nozzle, kappa, PbP0);
        let mut quans = calc_Mx(&nozzle, kappa, PbP0, nagare);

        //結果の出力
        let title = format!("./out/{}.png", i);
        let root = BitMapBackend::new(&title, (800, 600)).into_drawing_area();
        root.fill(&WHITE)?;
        let mut chart = ChartBuilder::on(&root)
            .margin(20)
            .x_label_area_size(10)
            .y_label_area_size(10)
            .build_cartesian_2d(-1f64..50f64, -1.0..50f64)?;
        
        chart
            .configure_mesh()
            .disable_x_mesh()
            .disable_y_mesh()
            .draw()?;
        let plotarea = chart.plotting_area();

        for q in quans {
            let mut r = 0.0;
            while r < ( q.rx * 1000.0 ) {
                plotarea.draw_pixel((q.x * 1000.0, r), &HSLColor( 0.50 +0.1*q.Mx, 1.0, 0.45));
                r += 0.01;
            }
        }

        P0 += 1.0e3;
        i += 1;
    }
    */

    Ok(())
}
