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
    println!("Me = {}", M);
    return M;
}
fn calc_M_super_from_AAc(kappa : &f64, AAc:&f64) -> f64{
    let mut M :f64 = 1.0;
    while *AAc - ( 1.0 / M * (( (kappa - 1.0) * M.powf(2.0) + 2.0) / (kappa + 1.0)).powf( (kappa + 1.0 )/( 2.0*(kappa - 1.0) )))
     > 1.0e-6
    {
        M += 0.01;
    }
    println!("超音速のM = {:.3}", M);
    return M;
}
fn calc_M_sub_from_AAc(kappa : &f64, AAc:&f64) -> f64{
    let mut M :f64 = 0.0;
    while *AAc - ( 1.0 / M * (( (kappa - 1.0) * M.powf(2.0) + 2.0) / (kappa + 1.0)).powf( (kappa + 1.0 )/( 2.0*(kappa - 1.0) )))
     > 1.0e-6
    {
        M += 0.01;
    }
    println!("亜音速のM = {:.3}", M);
    return M;
}
fn main() {
    let mut nozzle = Nozzle::new(1.0e-2, 1.5e-2, 2.0e-2);
    nozzle.init();
    //let AeAc = &nozzle.Ae / &nozzle.Ac;
    let AeAc = 2.403;

    let kappa : f64 = 1.40;
    let Mw = 28.8;
    let RR = 8.3143e3;//J/kmolK
    let R = RR/Mw;
    println!("R = {:.3}[kJ/kgK]", R/1000.0);

    let P0 : f64 = 600.0e3;//Pa
    let T0 : f64 = 290.0;// + 273.15;//K
    let Pa = 102.0e3;//Pa

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

    if pdp0 < PbP0 && PbP0 < 1.0 {
        pep0 = PbP0;
        Me = calc_M_sub_from_AAc(&kappa, &AeAc);
        println!("全域で亜音速");
    }
    else if PbP0 == pdp0 {
        pep0 = PbP0;
        Me = calc_M_sub_from_AAc(&kappa, &AeAc);
        println!("スロートで音速、その他で亜音速");
    }
    else if php0 < PbP0 && PbP0 < pdp0 {
        println!("非等エントロピー");
        println!("ノズルの末広部に垂直衝撃波");
    }
    else if php0 == PbP0 {
        Me = calc_M_super_from_AAc(&kappa, &AeAc);
        println!("出口に垂直衝撃波");
    }
    else if pjp0 < PbP0 && PbP0 < php0 {
        Me = calc_M_super_from_AAc(&kappa, &AeAc);
        println!("過膨張噴流");
    }
    else if pjp0 == PbP0 {
        Me = calc_M_super_from_AAc(&kappa, &AeAc);
        println!("適正膨張");
    } 
    else if PbP0 < pjp0 {
        Me = calc_M_super_from_AAc(&kappa, &AeAc);
        println!("不足膨張");
    }

    println!("pe = {:.2}[kPa]", pe/1000.0);
    let ue = ( 2.0 * kappa / ( kappa - 1.0) * R * T0 * ( 1.0 - ( pe / P0 ).powf( ( kappa - 1.0 ) / kappa ) ) ).powf(0.5);
    let tc = Tc * ( kappa + 1.0 ) / 2.0;

    let ae = ue / Me;

    println!("ue = {:.2}[m/s]     ae = {:.2}[m/s]    Me = {:.2}", ue, ae, Me);
}
