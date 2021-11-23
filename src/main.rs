use std::fs::File;
use std::io::{self, BufRead, Write, BufReader};

struct Nozzle{
    Dc : f64,
    De : f64,
    Le : f64,
    Ae : f64,
    Ac : f64,
}

impl Nozzle{
    fn calcAx(&self, x: f64) -> f64{
        let Dx = (self.De - self.Dc) * x / self.Le + self.Dc;
        return std::f64::consts::PI / 4.0 * Dx.powf(2.0); 
    }
    fn new_by_AeAc(AeAc:f64) -> Nozzle {
        Nozzle {Dc : 0.0, De: 0.0, Le: 0.0, Ae: AeAc, Ac: 1.0}
    }
    fn new(Dcm:f64, Dem:f64, Lem:f64) -> Nozzle {
        Nozzle {Dc : Dcm, De: Dem, Le: Lem, Ae: 0.0, Ac: 0.0}
    }
    fn calc_x (&self, AAc: &f64) -> f64 {
        let mut x = 0.0;
        while ( self.calcAx(x) / self.Ac - AAc ) < 1.0e-6
        {
            x += 0.001;
        }  
        println!("x = {:.3}[m]", x);
        return x;
    }

    fn init(&mut self){
        self.Ae = std::f64::consts::PI / 4.0 * self.De.powf(2.0); 
        self.Ac = std::f64::consts::PI / 4.0 * self.Dc.powf(2.0); 
        println!("Ae = {:.2}[mm2]   A* = {:.2}[mm2]", self.Ae*1.0e6, self.Ac*1.0e6);
    }
}

fn calcPc(P0: &f64, gammma: &f64) -> f64 {
    let pc = P0*( 2.0 / (gammma + 1.0) ).powf(gammma / (gammma - 1.0));
    println!("p* = {:.2}[kPa]", pc/1000.0);
    return pc;
}
fn calcTc(T0: &f64, gammma : &f64 ) -> f64{
    let Tc = T0*( 2.0 / (gammma + 1.0) );
    println!("T* = {:.2}[K]", Tc);
    return Tc;
}

fn calc_pd_by_p0(kappa: &f64, AeAc: &f64) -> f64{
    let mut pdp0:f64 = 0.528;
    while AeAc - ( (kappa - 1.0) / 2.0 * ( 2.0 / ( kappa + 1.0 ) ).powf( ( kappa + 1.0 )/( kappa - 1.0) )
            / ( pdp0.powf(2.0 / kappa ) - pdp0.powf( ( kappa + 1.0 )/kappa ) )  ).powf(0.5)
            > 1.0e-6
    {
        pdp0 += 0.01;
    }
    println!("pd/p0 = {:.3}", pdp0);
    return pdp0;
}
fn calc_pj_by_p0(kappa: &f64, AeAc: &f64) -> f64{
    let mut pjp0:f64 = 0.528;
    while AeAc - ( (kappa - 1.0) / 2.0 * ( 2.0 / ( kappa + 1.0 ) ).powf( ( kappa + 1.0 )/( kappa - 1.0) )
            / ( pjp0.powf(2.0 / kappa ) - pjp0.powf( ( kappa + 1.0 )/kappa ) )  ).powf(0.5)
             > 1.0e-6
    {
        pjp0 -= 0.01;
    }
    println!("pj/p0 = {:.3}", pjp0);
    return pjp0;
}
fn calc_ph_by_p0 (kappa: &f64, PbP0: &f64) -> f64{
    let Me = (2.0/(kappa - 1.0)*(PbP0.powf((1.0 - kappa)/kappa) -1.0)).powf(0.5);

    let php0 = (2.0 * kappa * Me.powf(2.0) - (kappa - 1.0))/(kappa + 1.0)*(1.0 + (kappa - 1.0)/2.0 * Me.powf(2.0)).powf(-(kappa/ (kappa - 1.0)));
    println!("ph/p0 = {:.3}", php0);
    return php0;
}
fn calc_M (kappa : &f64, pp0: &f64) -> f64 {
    let M = (2.0/(kappa - 1.0)*(pp0.powf((1.0 - kappa)/kappa) -1.0)).powf(0.5);
    println!("Me = {:.3}", M);
    return M;
} 
fn calc_M_before_shock(kappa: &f64, PeP0: &f64) -> f64 {
    let mut M1 : f64 = 1.0;
    while (
            ((kappa + 1.0) * M1.powf(2.0) / ((kappa - 1.0) * M1.powf(2.0) + 2.0)).powf(kappa / (kappa - 1.0))
          * ((kappa + 1.0)/(2.0 * kappa *  M1.powf(2.0) - (kappa - 1.0) )).powf( 1.0/(kappa -1.0) )
          ) - PeP0 > 1.0e-6
    {
        M1 += 0.01;   
    }
    println!("上流M1 = {:.3}", M1);
    return M1;
}
fn calc_M_after_shock(kappa: &f64, M1: &f64) -> f64 {
    let M2 = ( ((kappa -1.0) * M1.powf(2.0) + 2.0) / (2.0 * kappa * M1.powf(2.0) - (kappa - 1.0)) ).powf(0.5);
    println!("下流M2 = {:.3}", M2);
    return M2;   
}
fn calc_sub_sonic_M2_from_M1_by_AA(kappa: &f64, M1: &f64, A2A1: &f64) -> f64 {
    let mut M2 :f64 = 1.0;
    while (
         M1 / M2 * ( ( (kappa - 1.0) * M2.powf(2.0) + 2.0 )/( (kappa - 1.0) * M1.powf(2.0) + 2.0 ) ).powf( (kappa + 1.0)/(2.0 * (kappa - 1.0)) )
         ) - A2A1 < 1.0e-6
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
         ) - A2A1 < 1.0e-6
    {
        M2 += 0.001;
    }
    //println!("M2 = {:.3}", M2);
    return M2;
}
fn calc_AAc(kappa : &f64, M: &f64) -> f64 {
    let AAc = 1.0 / M * ( ((kappa - 1.0)*M.powf(2.0) + 2.0)/ (kappa + 1.0)).powf( (kappa + 1.0)/(2.0 * (kappa - 1.0)) );
    println!("A/Ac = {:.3}", AAc);
    return AAc;
}

fn main() -> Result<(), Box<std::error::Error>>{
    let mut nozzle = Nozzle::new(1.0e-2, 1.5e-2, 2.0e-2);
    nozzle.init();
    let AeAc = &nozzle.Ae / &nozzle.Ac;
    //let AeAc = 2.403;

    let kappa : f64 = 1.40;
    let Mw = 28.8;
    let RR = 8.3143e3;//J/kmolK
    let R = RR/Mw;
    println!("R = {:.3}[kJ/kgK]", R/1000.0);

    let P0 : f64 = 160.0e3;//Pa
    let T0 : f64 = 290.0;// + 273.15;//K
    let Pa = 101.325e3;//Pa

    let pc = calcPc(&P0, &kappa);
    let Tc = calcTc(&T0, &kappa);
    let Rowc = P0 / (R*T0);

    let PbP0 = Pa / P0;
    let mut pep0 = 0.0;
    let mut pe = 0.0;

    let pdp0 = calc_pd_by_p0(&kappa, &AeAc); 
    let pjp0 = calc_pj_by_p0(&kappa, &AeAc); 
    let php0 = calc_ph_by_p0(&kappa, &PbP0);      

    let mut Me = 0.0;
    let mut Mth = 0.0;

    let mut isIsentropic = true;

    //衝撃波がある場合
    let mut AmAc = 0.0; 

    if pdp0 <= PbP0 && PbP0 < 1.0 {
        //亜音速流
        pep0 = PbP0;
        if (pdp0 < PbP0 && PbP0 < 1.0){
            println!("全域で亜音速");
            //******************************
            //FIXME             Mth = 0.1;
            //*****************************
        }
        else {
            println!("スロートで音速、その他で亜音速");
            Mth = 1.0;
        }
        Me = calc_sub_sonic_M2_from_M1_by_AA(&kappa, &Mth, &AeAc);
    }
    else{ 
        //スロートで音速
        Mth = 1.0;
        Me = calc_super_sonic_M2_from_M1_by_AA(&kappa, &Mth, &AeAc);//衝撃波がある場合は上書きされる
        if php0 < PbP0 && PbP0 < pdp0 {
            println!("非等エントロピー");
            println!("ノズルの末広部に垂直衝撃波");
            isIsentropic = false;

            let M1 = calc_M_before_shock(&kappa, &PbP0);
            let M2 = calc_M_after_shock(&kappa, &M1);
            AmAc = calc_AAc(&kappa, &M1);//衝撃波発生位置の断面積比
            let xM = nozzle.calc_x(&AmAc);//衝撃波の発生位置
            let AeAm = AeAc / AmAc;
            Me = calc_sub_sonic_M2_from_M1_by_AA(&kappa, &M2, &AeAm);
        }
        else if php0 == PbP0 { 
            println!("出口に垂直衝撃波");
        }
        else if pjp0 < PbP0 && PbP0 < php0 {
            println!("過膨張噴流");
        }
        else if pjp0 == PbP0 {
            println!("適正膨張");
        } 
        else if PbP0 < pjp0 {
            println!("不足膨張");
        }
    }

    //ここから格子計算
    const N : usize =  256 * 2;
    let mut Mac_Array : [f64; N] = [0.0; N];
    Mac_Array[0] = Mth;
    let dx = &nozzle.Le / N as f64;

    for i in 1..N {
        let AxAc = nozzle.calcAx(dx * i as f64) / nozzle.Ac;
        let mut kyokushoM = 0.0; 
        if Me < 1.0 && isIsentropic{
            //等エントロピー、出口で亜音速
            kyokushoM = calc_sub_sonic_M2_from_M1_by_AA(&kappa, &Mth, &AxAc);
        }
        else if isIsentropic {
            //等エントロピー、出口で超音速
            kyokushoM = calc_super_sonic_M2_from_M1_by_AA(&kappa, &Mth, &AxAc);
        }
        else {
            //非等エントロピー
            if AxAc < AmAc {
                //衝撃波の上流 
                kyokushoM = calc_super_sonic_M2_from_M1_by_AA(&kappa, &Mth, &AxAc);
            }
            else {
                //衝撃波の下流
                let AxAm = AxAc / AmAc;
                kyokushoM = calc_sub_sonic_M2_from_M1_by_AA(&kappa, &Mth, &AxAm);            
            }
        }

        Mac_Array[i] = kyokushoM;
        //print!("{},", A2A1);
        //println!("{}",i);
    }

    //結果の出力
    println!("Me = {:.3}", Me);
    let mut f = std::fs::File::create(format!("./{:.3}_{:.3}_{:.3}.csv", Mth, Me, N))?;
    for i in 0 .. N {
       writeln!(f, "{},{}",i, Mac_Array[i])?;
    }
    f.flush()?;
    Ok(())
}
