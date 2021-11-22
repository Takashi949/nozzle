fn calcPc(P0: &f64, gammma: &f64) -> f64 {
    let pc = P0*( 2.0 / (gammma + 1.0) ).powf(( gammma / (gammma - 1.0) ));
    println!("p* = {:.2}[kPa]", pc/1000.0);
    return pc;
}
fn calcTc(T0: &f64, gammma : &f64 ) -> f64{
    let Tc = T0*( 2.0 / (gammma + 1.0) );
    println!("T* = {:.2}[K]", Tc);
    return Tc;
}

fn main() {
    let mut inpt = String::new();
    std::io::stdin().read_line(&mut inpt);
    let nozzleType = (inpt.starts_with("y"));//true => con-di nozzle
 
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
    let Rowc = P0/(R*T0);

    let mut pe = 0.0;
    if (Pa > pc) {
        pe = Pa;
        println!("亜音速噴流");
    }
    else if (Pa == pc){
        pe = pc;
        println!("音速噴流");
    }
    else if (Pa < pc){
        pe = pc;
        println!("不足膨張噴流");
    }
    println!("pe = {:.2}[kPa]", pe/1000.0);
    let ue = ( 2.0 * kappa / ( kappa - 1.0) * R * T0 * ( 1.0 - ( pe / P0 ).powf( ( kappa - 1.0 ) / kappa ) ) ).powf(0.5);
    let tc = Tc * ( kappa + 1.0 ) / 2.0;

    let Me = (2.0 / ( kappa - 1.0 ) *( ( P0 / pe ).powf( (kappa - 1.0 ) / kappa ) - 1.0 ) ).powf(0.5);
    let ae = ue / Me;

    println!("ue = {:.2}[m/s]     ae = {:.2}[m/s]    Me = {:.2}", ue, ae, Me);
}
