int flip_B3_B3(int phase){
  switch(phase){
  case 1:
  case 4:
    return phase; //cannot flip CC/RR
  case 2:
    return 3;
  case 3:
    return 2;
  default:
    return -1; //error
  }
  return -1; //error
}

int flip_B3_D1(int phase){
  switch(phase){
  case 1:
    return 2;
  case 2:
    return 1;
  case 3:
    return 4;
  case 4:
    return 3;
  default:
    return -1; //error
  }
  return -1; //error
}

int flip_B3_D2(int phase){
  switch(phase){
  case 1:
    return 3;
  case 2:
    return 4;
  case 3:
    return 1;
  case 4:
    return 2;
  default:
    return -1; //error
  }
  return -1; //error
}

int flip_D1_D1(int phase){
  switch(phase){
  case 1:
    return 2;
  case 2:
    return 1;
  case 3:
    return 4;
  case 4:
    return 3;
  default:
    return -1; //error
  }
  return -1; //error
}

int flip_D2_D2(int phase){
  switch(phase){
  case 1:
    return 3;
  case 2:
    return 4;
  case 3:
    return 1;
  case 4:
    return 2;
  default:
    return -1; //error
  }
  return -1; //error
}
