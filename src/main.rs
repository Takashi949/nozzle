use std::{io::Write, num::NonZeroI128, ops::Add};
use plotters::prelude::*;

const ERROR :f64 = 1.0e-3;
const delta :f64 = 1.0e-4;
const N:u64 = 1024;

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
        Nozzle {Di, Dth, xth, De, xe, ..Default::default()}
    }
    fn new_by_AeAc(AeAc:f64) -> Nozzle {
        Nozzle {Ae: AeAc, Ath: 1.0, ..Default::default()}
    }

    //先細部
    fn calc_AxC(&self, x: &f64) -> f64 {
        return std::f64::consts::PI / 4.0 * (self.calc_DxC(x)).powf(2.0); 
    }
    fn calc_DxC(&self, x:&f64) -> f64{
        let Dx = (self.Dth - self.Di) * x / self.xth + self.Di;
        return Dx;
    }
    //末広部
    fn calc_AxD(&self, x: &f64) -> f64{
        return std::f64::consts::PI / 4.0 * self.calc_DxD(&x).powf(2.0); 
    }
    fn calc_DxD(&self, x: &f64) -> f64{
        let x_xth = x - self.xth; 
        let Dx = (self.De - self.Dth) * x_xth / self.xe + self.Dth;
        return Dx; 
    }

    fn calc_xD (&self, AAc: &f64) -> f64 {
        let mut x = 0.0;
        while ( self.calc_AxD(&x) / self.Ath - AAc ) < ERROR
        {
            x += 0.001;
        }  
        println!("x = {:.3}[m]", x);
        return x;
    }

    fn init(&mut self){
        self.Ai = std::f64::consts::PI / 4.0 * self.Di.powf(2.0);
        self.Ae = std::f64::consts::PI / 4.0 * self.De.powf(2.0); 
        self.Ath = std::f64::consts::PI / 4.0 * self.Dth.powf(2.0); 
        println!("Ae = {:.2}[mm2]   Ath = {:.2}[mm2]", self.Ae*1.0e6, self.Ath*1.0e6);
    }
}
fn calc_sub_sonic_M2_from_M1_by_AA(kappa: &f64, M1: &f64, A2A1: &f64) -> f64 {
    let mut M2 :f64 = 1.0;
    while (
         M1 / M2 * ( ( (kappa - 1.0) * M2.powf(2.0) + 2.0 )/( (kappa - 1.0) * M1.powf(2.0) + 2.0 ) ).powf( (kappa + 1.0)/(2.0 * (kappa - 1.0)) )
         ) - A2A1 < ERROR
    {
        M2 -= 0.001;
    }
    //println!("M2 = {:.3}", M2);
    return M2;
}
fn calc_super_sonic_M2_from_M1_by_AA(kappa: &f64, M1: &f64, A2A1: &f64) -> f64 {
    let mut M2 :f64 = 1.0;
    while (
         M1 / M2 * ( ( (kappa - 1.0) * M2.powf(2.0) + 2.0 )/( (kappa - 1.0) * M1.powf(2.0) + 2.0 ) ).powf( (kappa + 1.0)/(2.0 * (kappa - 1.0)) )
         ) - A2A1 < ERROR
    {
        M2 += 0.001;
    }
    //println!("M2 = {:.3}", M2);
    return M2;
}
fn calc_M_before_shock(kappa: &f64, PeP0: &f64) -> f64 {
    //(8.18)
    let mut M1 : f64 = 1.0;
    while (
            ((kappa + 1.0) * M1.powf(2.0) / ((kappa - 1.0) * M1.powf(2.0) + 2.0)).powf(kappa / (kappa - 1.0))
          * ((kappa + 1.0)/(2.0 * kappa *  M1.powf(2.0) - (kappa - 1.0) )).powf( 1.0/(kappa -1.0) )
          ) - PeP0 > ERROR
    {
        M1 += 0.01;   
    }
    println!("上流M1 = {:.3}", M1);
    return M1;
}
fn calc_M_after_shock(kappa: &f64, M1: &f64) -> f64 {
    //(8.18)
    let M2 = ( ((kappa -1.0) * M1.powf(2.0) + 2.0) / (2.0 * kappa * M1.powf(2.0) - (kappa - 1.0)) ).powf(0.5);
    println!("下流M2 = {:.3}", M2);
    return M2;   
}
fn calc_AAc(kappa : &f64, M: &f64) -> f64 {
    let AAc = 1.0 / M * ( ((kappa - 1.0)*M.powf(2.0) + 2.0)/ (kappa + 1.0)).powf( (kappa + 1.0)/(2.0 * (kappa - 1.0)) );
    println!("A/Ac = {:.3}", AAc);
    return AAc;
}
fn siki_7_21(kappa : &f64, pdp0: &f64) -> f64{
    //pdp0からAe/A*=
    ( ( (kappa - 1.0) / 2.0 )*( 2.0 / ( kappa + 1.0 ) ).powf( ( kappa + 1.0 )/( kappa - 1.0) )
    / ( pdp0.powf(2.0 / kappa ) - pdp0.powf( ( kappa + 1.0 )/kappa ) )  ).powf(0.5)
}
fn siki_8_25(kappa : &f64, M1: f64) -> f64 {
    //p2: 衝撃は直後の静圧 p0:ノズル上流全圧
    let p2p0 = (2.0 * kappa * M1.powf(2.0) - (kappa - 1.0))/(kappa + 1.0)*(1.0 + (kappa - 1.0)/2.0 * M1.powf(2.0)).powf(-(kappa/ (kappa - 1.0)));
    return p2p0;
}

fn calc_pcp0(kappa : &f64) -> f64 {
    let pcp0 = (2.0 / (kappa + 1.0)).powf(kappa / (kappa - 1.0));
    println!("p*/p0 = {:.3}", pcp0);
    return pcp0;
}

fn calcPc(P0: &f64, gammma: &f64) -> f64 {
    let pc = P0* calc_pcp0(gammma);
    println!("p* = {:.2}[kPa]", pc/1000.0);
    return pc;
}
fn calcTc(T0: &f64, gammma : &f64 ) -> f64{
    let Tc = T0*( 2.0 / (gammma + 1.0) );
    println!("T* = {:.2}[K]", Tc);
    return Tc;
}

fn calc_pd_by_p0(kappa: &f64, AeAc: &f64) -> f64{
    let mut pdp0:f64 = calc_pcp0(&kappa);
    while ( AeAc - siki_7_21(kappa, &pdp0) ).abs() > ERROR
    {
        pdp0 += delta;
    }
    println!("pd/p0 = {:.3}", pdp0);
    return pdp0;
}
fn calc_pj_by_p0(kappa: &f64, AeAc: &f64) -> f64{
    let mut pjp0:f64 = calc_pcp0(&kappa);
    while (AeAc - siki_7_21(kappa, &pjp0) ).abs() > ERROR
    {
        pjp0 -= delta;
    }
    println!("pj/p0 = {:.3}", pjp0);
    return pjp0;
}
fn calc_ph_by_p0 (kappa: &f64, AeAc: &f64) -> f64{
    let Me = 1.0;
    let mut php0 = siki_8_25(&kappa, calc_super_sonic_M2_from_M1_by_AA(&kappa, &Me, AeAc)); 
    println!("ph/p0 = {:.3}", php0);
    return php0;
}
fn calc_M (kappa : &f64, pp0: &f64) -> f64 {
    let M = (2.0/(kappa - 1.0)*(pp0.powf((1.0 - kappa)/kappa) -1.0)).powf(0.5);
    println!("M = {:.3}", M);
    return M;
} 

fn check_nagare_keitai(nozzle : &Nozzle, kappa : &f64, PbP0: f64) -> NAGARE{
    let AeAth = &nozzle.Ae / &nozzle.Ath;
    let pdp0 = calc_pd_by_p0(&kappa, &AeAth); 
    let pjp0 = calc_pj_by_p0(&kappa, &AeAth); 
    let php0 = calc_ph_by_p0(&kappa, &AeAth);   

    if pdp0 <= PbP0 && PbP0 < 1.0 {
        //亜音速流
        if (pdp0 < PbP0 && PbP0 < 1.0){
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
        if php0 < PbP0 && PbP0 < pdp0 {
            println!("非等エントロピー");
            println!("ノズルの末広部に垂直衝撃波");
            return NAGARE::suehiro;
        }
        else if php0 == PbP0 { 
            println!("出口に垂直衝撃波");
            return NAGARE::deguti;
        }
        else if pjp0 < PbP0 && PbP0 < php0 {
            println!("過膨張噴流");
            return NAGARE::ka;
        }
        else if pjp0 == PbP0 {
            println!("適正膨張");
            return NAGARE::saiteki;
        } 
        else if PbP0 < pjp0 {
            println!("不足膨張");
            return NAGARE::husoku;       
        }
    }
    return NAGARE::humei;
}

fn calc_Mx(nozzle : &Nozzle, kappa : &f64, PbP0: f64, nagare: NAGARE) -> Vec<quantity>{
    let AeAth = &nozzle.Ae / &nozzle.Ath;
    let mut Mth = 1.0;//MthL1のときだけ変更
    let mut Me = 0.0;

    //衝撃波がある場合
    let mut AmAc= 0.0;

    match nagare {
        NAGARE::MthL1 => {
            Mth = calc_M(&kappa, &PbP0);
            Me = calc_sub_sonic_M2_from_M1_by_AA(&kappa, &Mth, &AeAth);
        },
        NAGARE::MthG1MeL1 => {
            Me = calc_sub_sonic_M2_from_M1_by_AA(&kappa, &Mth, &AeAth);
        },
        NAGARE::suehiro => {
            let M1 = calc_M_before_shock(&kappa, &PbP0);
            let M2 = calc_M_after_shock(&kappa, &M1);
            AmAc = calc_AAc(&kappa, &M1);//衝撃波発生位置の断面積比
            let xM = nozzle.calc_xD(&AmAc);//衝撃波の発生位置
            let AeAm = AeAth / AmAc;
            Me = calc_sub_sonic_M2_from_M1_by_AA(&kappa, &M2, &AeAm);
        },
        _=>{
            Me = calc_super_sonic_M2_from_M1_by_AA(&kappa, &Mth, &AeAth);//衝撃波がある場合は上書きされる
        }
    }
    println!("Me = {:.3}    Mth = {:.3}",Me, Mth);

    //ここから格子計算
    let mut Quans : Vec<quantity> = Vec::new();
    let dx = &nozzle.xe / N as f64;

    for i in 1..N {
        let x = dx * i as f64;
        let mut Mx = 0.0;
        let mut rx = 0.0;
        if x <= nozzle.xth {
            //先細部の計算
            // 状態2:x     状態1:th
            let AxAth = nozzle.calc_AxC(&x) / nozzle.Ath;
            rx = nozzle.calc_DxC(&x) / 2.0;
            Mx = calc_sub_sonic_M2_from_M1_by_AA(&kappa, &Mth, &AxAth);
        }
        else {
            //末広部の計算
            let AxAth = nozzle.calc_AxD(&x) / nozzle.Ath;
            rx = nozzle.calc_DxD(&x) / 2.0;

            if let NAGARE::suehiro = nagare{
                //非等エントロピー
                if AxAth < AmAc {
                    //衝撃波の上流 
                    Mx = calc_super_sonic_M2_from_M1_by_AA(&kappa, &Mth, &AxAth);
                }
                else {
                    //衝撃波の下流
                    let AxAm = AxAth / AmAc;
                    Mx = calc_sub_sonic_M2_from_M1_by_AA(&kappa, &Mth, &AxAm);            
                }
            }
            else if Me < 1.0{
                //等エントロピー、亜音速
                Mx = calc_sub_sonic_M2_from_M1_by_AA(&kappa, &Mth, &AxAth);
            }
            else {
                //等エントロピー、超音速
                Mx = calc_super_sonic_M2_from_M1_by_AA(&kappa, &Mth, &AxAth);
            }
        }

        //Mac_Array[i] = kyokushoM;
        //print!("{},", A2A1);
        //println!("{}",i);

        let res = quantity {x, Mx, rx};
        Quans.push(res);
    }
    return Quans;
}
fn main() -> Result<(), Box<std::error::Error>>{
    let mut nozzle = Nozzle::new(20.0e-3, 10.0e-3, 10.0e-3, 17.320508e-3, 40.0e-3);
    nozzle.init();
    let AeAth = &nozzle.Ae / &nozzle.Ath;
    //let AeAc = 2.403;
    println!("AeAth = {:.3}",AeAth);

    let kappa : f64 = 1.40;
    let Mw = 28.8;
    let RR = 8.3143e3;//J/kmolK
    let R = RR/Mw;
    println!("R = {:.3}[kJ/kgK]", R/1000.0);

    //let P0 : f64 = 102.0e3;//全域亜音速
    let P0 = 150.0e3;//衝撃波
    let T0 : f64 = 290.0;// + 273.15;//K
    let Pa = 101.325e3;//Pa

    let PbP0 = Pa / P0;
    println!("PbP0 = {:.3}", PbP0);

    let mut Me = 0.0;
    let mut Mth = 0.0;

    let nagare = check_nagare_keitai(&nozzle, &kappa, PbP0);
    let mut quans = calc_Mx(&nozzle, &kappa, PbP0, nagare);

    //結果の出力
    let root = BitMapBackend::new("./out.png", (800, 600)).into_drawing_area();
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
    let range = plotarea.get_pixel_range();

    for q in quans {
        let mut r = 0.0;
        while r < ( q.rx * 1000.0 ) {
            plotarea.draw_pixel((q.x * 1000.0, r), &HSLColor( 0.50 +0.3*q.Mx, 1.0, 0.5));
            r += 0.001;
        }
    }
    /*
    let mut f = std::fs::File::create(format!("./{}.csv", PbP0))?;
    for q in quans {
       writeln!(f, "{},{}, {}", q.x, q.Mx, q.rx)?;
    }
    f.flush()?;
    */
    Ok(())
}
