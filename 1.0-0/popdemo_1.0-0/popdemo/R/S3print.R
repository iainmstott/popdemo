print.projection<-function(x,...){
    attributes(x)$class<-NULL
    print(x)
}

print.tfa<-function(x,...){
    attributes(x)$class<-NULL
    print(x)
}

print.tfam<-function(x,...){
    attributes(x)$class<-NULL
    l<-length(x)
    print(x[1:(l-1)])
}
