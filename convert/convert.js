null;var Module=typeof Module!=="undefined"?Module:{};var moduleOverrides={};var key;for(key in Module){if(Module.hasOwnProperty(key)){moduleOverrides[key]=Module[key]}}var arguments_=[];var thisProgram="./this.program";var quit_=function(status,toThrow){throw toThrow};var ENVIRONMENT_IS_WEB=typeof window==="object";var ENVIRONMENT_IS_WORKER=typeof importScripts==="function";var ENVIRONMENT_IS_NODE=typeof process==="object"&&typeof process.versions==="object"&&typeof process.versions.node==="string";var scriptDirectory="";function locateFile(path){if(Module["locateFile"]){return Module["locateFile"](path,scriptDirectory)}return scriptDirectory+path}var read_,readAsync,readBinary,setWindowTitle;var nodeFS;var nodePath;if(ENVIRONMENT_IS_NODE){if(ENVIRONMENT_IS_WORKER){scriptDirectory=require("path").dirname(scriptDirectory)+"/"}else{scriptDirectory=__dirname+"/"}read_=function shell_read(filename,binary){if(!nodeFS)nodeFS=require("fs");if(!nodePath)nodePath=require("path");filename=nodePath["normalize"](filename);return nodeFS["readFileSync"](filename,binary?null:"utf8")};readBinary=function readBinary(filename){var ret=read_(filename,true);if(!ret.buffer){ret=new Uint8Array(ret)}assert(ret.buffer);return ret};readAsync=function readAsync(filename,onload,onerror){if(!nodeFS)nodeFS=require("fs");if(!nodePath)nodePath=require("path");filename=nodePath["normalize"](filename);nodeFS["readFile"](filename,function(err,data){if(err)onerror(err);else onload(data.buffer)})};if(process["argv"].length>1){thisProgram=process["argv"][1].replace(/\\/g,"/")}arguments_=process["argv"].slice(2);if(typeof module!=="undefined"){module["exports"]=Module}process["on"]("uncaughtException",function(ex){if(!(ex instanceof ExitStatus)){throw ex}});process["on"]("unhandledRejection",abort);quit_=function(status,toThrow){if(keepRuntimeAlive()){process["exitCode"]=status;throw toThrow}process["exit"](status)};Module["inspect"]=function(){return"[Emscripten Module object]"}}else if(ENVIRONMENT_IS_WEB||ENVIRONMENT_IS_WORKER){if(ENVIRONMENT_IS_WORKER){scriptDirectory=self.location.href}else if(typeof document!=="undefined"&&document.currentScript){scriptDirectory=document.currentScript.src}if(scriptDirectory.indexOf("blob:")!==0){scriptDirectory=scriptDirectory.substr(0,scriptDirectory.lastIndexOf("/")+1)}else{scriptDirectory=""}{read_=function(url){var xhr=new XMLHttpRequest;xhr.open("GET",url,false);xhr.send(null);return xhr.responseText};if(ENVIRONMENT_IS_WORKER){readBinary=function(url){var xhr=new XMLHttpRequest;xhr.open("GET",url,false);xhr.responseType="arraybuffer";xhr.send(null);return new Uint8Array(xhr.response)}}readAsync=function(url,onload,onerror){var xhr=new XMLHttpRequest;xhr.open("GET",url,true);xhr.responseType="arraybuffer";xhr.onload=function(){if(xhr.status==200||xhr.status==0&&xhr.response){onload(xhr.response);return}onerror()};xhr.onerror=onerror;xhr.send(null)}}setWindowTitle=function(title){document.title=title}}else{}var out=Module["print"]||console.log.bind(console);var err=Module["printErr"]||console.warn.bind(console);for(key in moduleOverrides){if(moduleOverrides.hasOwnProperty(key)){Module[key]=moduleOverrides[key]}}moduleOverrides=null;if(Module["arguments"])arguments_=Module["arguments"];if(Module["thisProgram"])thisProgram=Module["thisProgram"];if(Module["quit"])quit_=Module["quit"];var tempRet0=0;var setTempRet0=function(value){tempRet0=value};var getTempRet0=function(){return tempRet0};var wasmBinary;if(Module["wasmBinary"])wasmBinary=Module["wasmBinary"];var noExitRuntime=Module["noExitRuntime"]||true;if(typeof WebAssembly!=="object"){abort("no native wasm support detected")}function getValue(ptr,type,noSafe){type=type||"i8";if(type.charAt(type.length-1)==="*")type="i32";switch(type){case"i1":return HEAP8[ptr>>0];case"i8":return HEAP8[ptr>>0];case"i16":return HEAP16[ptr>>1];case"i32":return HEAP32[ptr>>2];case"i64":return HEAP32[ptr>>2];case"float":return HEAPF32[ptr>>2];case"double":return HEAPF64[ptr>>3];default:abort("invalid type for getValue: "+type)}return null}var wasmMemory;var ABORT=false;var EXITSTATUS;function assert(condition,text){if(!condition){abort("Assertion failed: "+text)}}function getCFunc(ident){var func=Module["_"+ident];assert(func,"Cannot call unknown function "+ident+", make sure it is exported");return func}function ccall(ident,returnType,argTypes,args,opts){var toC={"string":function(str){var ret=0;if(str!==null&&str!==undefined&&str!==0){var len=(str.length<<2)+1;ret=stackAlloc(len);stringToUTF8(str,ret,len)}return ret},"array":function(arr){var ret=stackAlloc(arr.length);writeArrayToMemory(arr,ret);return ret}};function convertReturnValue(ret){if(returnType==="string")return UTF8ToString(ret);if(returnType==="boolean")return Boolean(ret);return ret}var func=getCFunc(ident);var cArgs=[];var stack=0;if(args){for(var i=0;i<args.length;i++){var converter=toC[argTypes[i]];if(converter){if(stack===0)stack=stackSave();cArgs[i]=converter(args[i])}else{cArgs[i]=args[i]}}}var ret=func.apply(null,cArgs);function onDone(ret){if(stack!==0)stackRestore(stack);return convertReturnValue(ret)}ret=onDone(ret);return ret}var UTF8Decoder=typeof TextDecoder!=="undefined"?new TextDecoder("utf8"):undefined;function UTF8ArrayToString(heap,idx,maxBytesToRead){var endIdx=idx+maxBytesToRead;var endPtr=idx;while(heap[endPtr]&&!(endPtr>=endIdx))++endPtr;if(endPtr-idx>16&&heap.subarray&&UTF8Decoder){return UTF8Decoder.decode(heap.subarray(idx,endPtr))}else{var str="";while(idx<endPtr){var u0=heap[idx++];if(!(u0&128)){str+=String.fromCharCode(u0);continue}var u1=heap[idx++]&63;if((u0&224)==192){str+=String.fromCharCode((u0&31)<<6|u1);continue}var u2=heap[idx++]&63;if((u0&240)==224){u0=(u0&15)<<12|u1<<6|u2}else{u0=(u0&7)<<18|u1<<12|u2<<6|heap[idx++]&63}if(u0<65536){str+=String.fromCharCode(u0)}else{var ch=u0-65536;str+=String.fromCharCode(55296|ch>>10,56320|ch&1023)}}}return str}function UTF8ToString(ptr,maxBytesToRead){return ptr?UTF8ArrayToString(HEAPU8,ptr,maxBytesToRead):""}function stringToUTF8Array(str,heap,outIdx,maxBytesToWrite){if(!(maxBytesToWrite>0))return 0;var startIdx=outIdx;var endIdx=outIdx+maxBytesToWrite-1;for(var i=0;i<str.length;++i){var u=str.charCodeAt(i);if(u>=55296&&u<=57343){var u1=str.charCodeAt(++i);u=65536+((u&1023)<<10)|u1&1023}if(u<=127){if(outIdx>=endIdx)break;heap[outIdx++]=u}else if(u<=2047){if(outIdx+1>=endIdx)break;heap[outIdx++]=192|u>>6;heap[outIdx++]=128|u&63}else if(u<=65535){if(outIdx+2>=endIdx)break;heap[outIdx++]=224|u>>12;heap[outIdx++]=128|u>>6&63;heap[outIdx++]=128|u&63}else{if(outIdx+3>=endIdx)break;heap[outIdx++]=240|u>>18;heap[outIdx++]=128|u>>12&63;heap[outIdx++]=128|u>>6&63;heap[outIdx++]=128|u&63}}heap[outIdx]=0;return outIdx-startIdx}function stringToUTF8(str,outPtr,maxBytesToWrite){return stringToUTF8Array(str,HEAPU8,outPtr,maxBytesToWrite)}function lengthBytesUTF8(str){var len=0;for(var i=0;i<str.length;++i){var u=str.charCodeAt(i);if(u>=55296&&u<=57343)u=65536+((u&1023)<<10)|str.charCodeAt(++i)&1023;if(u<=127)++len;else if(u<=2047)len+=2;else if(u<=65535)len+=3;else len+=4}return len}function writeArrayToMemory(array,buffer){HEAP8.set(array,buffer)}function writeAsciiToMemory(str,buffer,dontAddNull){for(var i=0;i<str.length;++i){HEAP8[buffer++>>0]=str.charCodeAt(i)}if(!dontAddNull)HEAP8[buffer>>0]=0}function alignUp(x,multiple){if(x%multiple>0){x+=multiple-x%multiple}return x}var buffer,HEAP8,HEAPU8,HEAP16,HEAPU16,HEAP32,HEAPU32,HEAPF32,HEAPF64;function updateGlobalBufferAndViews(buf){buffer=buf;Module["HEAP8"]=HEAP8=new Int8Array(buf);Module["HEAP16"]=HEAP16=new Int16Array(buf);Module["HEAP32"]=HEAP32=new Int32Array(buf);Module["HEAPU8"]=HEAPU8=new Uint8Array(buf);Module["HEAPU16"]=HEAPU16=new Uint16Array(buf);Module["HEAPU32"]=HEAPU32=new Uint32Array(buf);Module["HEAPF32"]=HEAPF32=new Float32Array(buf);Module["HEAPF64"]=HEAPF64=new Float64Array(buf)}var INITIAL_MEMORY=Module["INITIAL_MEMORY"]||16777216;var wasmTable;var __ATPRERUN__=[];var __ATINIT__=[];var __ATPOSTRUN__=[];var runtimeInitialized=false;var runtimeKeepaliveCounter=0;function keepRuntimeAlive(){return noExitRuntime||runtimeKeepaliveCounter>0}function preRun(){if(Module["preRun"]){if(typeof Module["preRun"]=="function")Module["preRun"]=[Module["preRun"]];while(Module["preRun"].length){addOnPreRun(Module["preRun"].shift())}}callRuntimeCallbacks(__ATPRERUN__)}function initRuntime(){runtimeInitialized=true;callRuntimeCallbacks(__ATINIT__)}function postRun(){if(Module["postRun"]){if(typeof Module["postRun"]=="function")Module["postRun"]=[Module["postRun"]];while(Module["postRun"].length){addOnPostRun(Module["postRun"].shift())}}callRuntimeCallbacks(__ATPOSTRUN__)}function addOnPreRun(cb){__ATPRERUN__.unshift(cb)}function addOnInit(cb){__ATINIT__.unshift(cb)}function addOnPostRun(cb){__ATPOSTRUN__.unshift(cb)}var runDependencies=0;var runDependencyWatcher=null;var dependenciesFulfilled=null;function addRunDependency(id){runDependencies++;if(Module["monitorRunDependencies"]){Module["monitorRunDependencies"](runDependencies)}}function removeRunDependency(id){runDependencies--;if(Module["monitorRunDependencies"]){Module["monitorRunDependencies"](runDependencies)}if(runDependencies==0){if(runDependencyWatcher!==null){clearInterval(runDependencyWatcher);runDependencyWatcher=null}if(dependenciesFulfilled){var callback=dependenciesFulfilled;dependenciesFulfilled=null;callback()}}}Module["preloadedImages"]={};Module["preloadedAudios"]={};function abort(what){{if(Module["onAbort"]){Module["onAbort"](what)}}what+="";err(what);ABORT=true;EXITSTATUS=1;what="abort("+what+"). Build with -s ASSERTIONS=1 for more info.";var e=new WebAssembly.RuntimeError(what);throw e}var dataURIPrefix="data:application/octet-stream;base64,";function isDataURI(filename){return filename.startsWith(dataURIPrefix)}function isFileURI(filename){return filename.startsWith("file://")}var wasmBinaryFile;wasmBinaryFile="convert.wasm";if(!isDataURI(wasmBinaryFile)){wasmBinaryFile=locateFile(wasmBinaryFile)}function getBinary(file){try{if(file==wasmBinaryFile&&wasmBinary){return new Uint8Array(wasmBinary)}if(readBinary){return readBinary(file)}else{throw"both async and sync fetching of the wasm failed"}}catch(err){abort(err)}}function getBinaryPromise(){if(!wasmBinary&&(ENVIRONMENT_IS_WEB||ENVIRONMENT_IS_WORKER)){if(typeof fetch==="function"&&!isFileURI(wasmBinaryFile)){return fetch(wasmBinaryFile,{credentials:"same-origin"}).then(function(response){if(!response["ok"]){throw"failed to load wasm binary file at '"+wasmBinaryFile+"'"}return response["arrayBuffer"]()}).catch(function(){return getBinary(wasmBinaryFile)})}else{if(readAsync){return new Promise(function(resolve,reject){readAsync(wasmBinaryFile,function(response){resolve(new Uint8Array(response))},reject)})}}}return Promise.resolve().then(function(){return getBinary(wasmBinaryFile)})}function createWasm(){var info={"a":asmLibraryArg};function receiveInstance(instance,module){var exports=instance.exports;Module["asm"]=exports;wasmMemory=Module["asm"]["U"];updateGlobalBufferAndViews(wasmMemory.buffer);wasmTable=Module["asm"]["la"];addOnInit(Module["asm"]["V"]);removeRunDependency("wasm-instantiate")}addRunDependency("wasm-instantiate");function receiveInstantiationResult(result){receiveInstance(result["instance"])}function instantiateArrayBuffer(receiver){return getBinaryPromise().then(function(binary){return WebAssembly.instantiate(binary,info)}).then(function(instance){return instance}).then(receiver,function(reason){err("failed to asynchronously prepare wasm: "+reason);abort(reason)})}function instantiateAsync(){if(!wasmBinary&&typeof WebAssembly.instantiateStreaming==="function"&&!isDataURI(wasmBinaryFile)&&!isFileURI(wasmBinaryFile)&&typeof fetch==="function"){return fetch(wasmBinaryFile,{credentials:"same-origin"}).then(function(response){var result=WebAssembly.instantiateStreaming(response,info);return result.then(receiveInstantiationResult,function(reason){err("wasm streaming compile failed: "+reason);err("falling back to ArrayBuffer instantiation");return instantiateArrayBuffer(receiveInstantiationResult)})})}else{return instantiateArrayBuffer(receiveInstantiationResult)}}if(Module["instantiateWasm"]){try{var exports=Module["instantiateWasm"](info,receiveInstance);return exports}catch(e){err("Module.instantiateWasm callback failed with error: "+e);return false}}instantiateAsync();return{}}function callRuntimeCallbacks(callbacks){while(callbacks.length>0){var callback=callbacks.shift();if(typeof callback=="function"){callback(Module);continue}var func=callback.func;if(typeof func==="number"){if(callback.arg===undefined){wasmTable.get(func)()}else{wasmTable.get(func)(callback.arg)}}else{func(callback.arg===undefined?null:callback.arg)}}}function ___assert_fail(condition,filename,line,func){abort("Assertion failed: "+UTF8ToString(condition)+", at: "+[filename?UTF8ToString(filename):"unknown filename",line,func?UTF8ToString(func):"unknown function"])}function ___cxa_allocate_exception(size){return _malloc(size+16)+16}function ExceptionInfo(excPtr){this.excPtr=excPtr;this.ptr=excPtr-16;this.set_type=function(type){HEAP32[this.ptr+4>>2]=type};this.get_type=function(){return HEAP32[this.ptr+4>>2]};this.set_destructor=function(destructor){HEAP32[this.ptr+8>>2]=destructor};this.get_destructor=function(){return HEAP32[this.ptr+8>>2]};this.set_refcount=function(refcount){HEAP32[this.ptr>>2]=refcount};this.set_caught=function(caught){caught=caught?1:0;HEAP8[this.ptr+12>>0]=caught};this.get_caught=function(){return HEAP8[this.ptr+12>>0]!=0};this.set_rethrown=function(rethrown){rethrown=rethrown?1:0;HEAP8[this.ptr+13>>0]=rethrown};this.get_rethrown=function(){return HEAP8[this.ptr+13>>0]!=0};this.init=function(type,destructor){this.set_type(type);this.set_destructor(destructor);this.set_refcount(0);this.set_caught(false);this.set_rethrown(false)};this.add_ref=function(){var value=HEAP32[this.ptr>>2];HEAP32[this.ptr>>2]=value+1};this.release_ref=function(){var prev=HEAP32[this.ptr>>2];HEAP32[this.ptr>>2]=prev-1;return prev===1}}function CatchInfo(ptr){this.free=function(){_free(this.ptr);this.ptr=0};this.set_base_ptr=function(basePtr){HEAP32[this.ptr>>2]=basePtr};this.get_base_ptr=function(){return HEAP32[this.ptr>>2]};this.set_adjusted_ptr=function(adjustedPtr){HEAP32[this.ptr+4>>2]=adjustedPtr};this.get_adjusted_ptr_addr=function(){return this.ptr+4};this.get_adjusted_ptr=function(){return HEAP32[this.ptr+4>>2]};this.get_exception_ptr=function(){var isPointer=___cxa_is_pointer_type(this.get_exception_info().get_type());if(isPointer){return HEAP32[this.get_base_ptr()>>2]}var adjusted=this.get_adjusted_ptr();if(adjusted!==0)return adjusted;return this.get_base_ptr()};this.get_exception_info=function(){return new ExceptionInfo(this.get_base_ptr())};if(ptr===undefined){this.ptr=_malloc(8);this.set_adjusted_ptr(0)}else{this.ptr=ptr}}var exceptionCaught=[];function exception_addRef(info){info.add_ref()}var uncaughtExceptionCount=0;function ___cxa_begin_catch(ptr){var catchInfo=new CatchInfo(ptr);var info=catchInfo.get_exception_info();if(!info.get_caught()){info.set_caught(true);uncaughtExceptionCount--}info.set_rethrown(false);exceptionCaught.push(catchInfo);exception_addRef(info);return catchInfo.get_exception_ptr()}var exceptionLast=0;function ___cxa_free_exception(ptr){return _free(new ExceptionInfo(ptr).ptr)}function exception_decRef(info){if(info.release_ref()&&!info.get_rethrown()){var destructor=info.get_destructor();if(destructor){wasmTable.get(destructor)(info.excPtr)}___cxa_free_exception(info.excPtr)}}function ___cxa_end_catch(){_setThrew(0);var catchInfo=exceptionCaught.pop();exception_decRef(catchInfo.get_exception_info());catchInfo.free();exceptionLast=0}function ___resumeException(catchInfoPtr){var catchInfo=new CatchInfo(catchInfoPtr);var ptr=catchInfo.get_base_ptr();if(!exceptionLast){exceptionLast=ptr}catchInfo.free();throw ptr}function ___cxa_find_matching_catch_2(){var thrown=exceptionLast;if(!thrown){setTempRet0(0);return 0|0}var info=new ExceptionInfo(thrown);var thrownType=info.get_type();var catchInfo=new CatchInfo;catchInfo.set_base_ptr(thrown);catchInfo.set_adjusted_ptr(thrown);if(!thrownType){setTempRet0(0);return catchInfo.ptr|0}var typeArray=Array.prototype.slice.call(arguments);for(var i=0;i<typeArray.length;i++){var caughtType=typeArray[i];if(caughtType===0||caughtType===thrownType){break}if(___cxa_can_catch(caughtType,thrownType,catchInfo.get_adjusted_ptr_addr())){setTempRet0(caughtType);return catchInfo.ptr|0}}setTempRet0(thrownType);return catchInfo.ptr|0}function ___cxa_find_matching_catch_3(){var thrown=exceptionLast;if(!thrown){setTempRet0(0);return 0|0}var info=new ExceptionInfo(thrown);var thrownType=info.get_type();var catchInfo=new CatchInfo;catchInfo.set_base_ptr(thrown);catchInfo.set_adjusted_ptr(thrown);if(!thrownType){setTempRet0(0);return catchInfo.ptr|0}var typeArray=Array.prototype.slice.call(arguments);for(var i=0;i<typeArray.length;i++){var caughtType=typeArray[i];if(caughtType===0||caughtType===thrownType){break}if(___cxa_can_catch(caughtType,thrownType,catchInfo.get_adjusted_ptr_addr())){setTempRet0(caughtType);return catchInfo.ptr|0}}setTempRet0(thrownType);return catchInfo.ptr|0}function ___cxa_find_matching_catch_4(){var thrown=exceptionLast;if(!thrown){setTempRet0(0);return 0|0}var info=new ExceptionInfo(thrown);var thrownType=info.get_type();var catchInfo=new CatchInfo;catchInfo.set_base_ptr(thrown);catchInfo.set_adjusted_ptr(thrown);if(!thrownType){setTempRet0(0);return catchInfo.ptr|0}var typeArray=Array.prototype.slice.call(arguments);for(var i=0;i<typeArray.length;i++){var caughtType=typeArray[i];if(caughtType===0||caughtType===thrownType){break}if(___cxa_can_catch(caughtType,thrownType,catchInfo.get_adjusted_ptr_addr())){setTempRet0(caughtType);return catchInfo.ptr|0}}setTempRet0(thrownType);return catchInfo.ptr|0}function ___cxa_rethrow(){var catchInfo=exceptionCaught.pop();if(!catchInfo){abort("no exception to throw")}var info=catchInfo.get_exception_info();var ptr=catchInfo.get_base_ptr();if(!info.get_rethrown()){exceptionCaught.push(catchInfo);info.set_rethrown(true);info.set_caught(false);uncaughtExceptionCount++}else{catchInfo.free()}exceptionLast=ptr;throw ptr}function ___cxa_throw(ptr,type,destructor){var info=new ExceptionInfo(ptr);info.init(type,destructor);exceptionLast=ptr;uncaughtExceptionCount++;throw ptr}function ___cxa_uncaught_exceptions(){return uncaughtExceptionCount}function _abort(){abort()}function _emscripten_memcpy_big(dest,src,num){HEAPU8.copyWithin(dest,src,src+num)}function emscripten_realloc_buffer(size){try{wasmMemory.grow(size-buffer.byteLength+65535>>>16);updateGlobalBufferAndViews(wasmMemory.buffer);return 1}catch(e){}}function _emscripten_resize_heap(requestedSize){var oldSize=HEAPU8.length;requestedSize=requestedSize>>>0;var maxHeapSize=2147483648;if(requestedSize>maxHeapSize){return false}for(var cutDown=1;cutDown<=4;cutDown*=2){var overGrownHeapSize=oldSize*(1+.2/cutDown);overGrownHeapSize=Math.min(overGrownHeapSize,requestedSize+100663296);var newSize=Math.min(maxHeapSize,alignUp(Math.max(requestedSize,overGrownHeapSize),65536));var replacement=emscripten_realloc_buffer(newSize);if(replacement){return true}}return false}var ENV={};function getExecutableName(){return thisProgram||"./this.program"}function getEnvStrings(){if(!getEnvStrings.strings){var lang=(typeof navigator==="object"&&navigator.languages&&navigator.languages[0]||"C").replace("-","_")+".UTF-8";var env={"USER":"web_user","LOGNAME":"web_user","PATH":"/","PWD":"/","HOME":"/home/web_user","LANG":lang,"_":getExecutableName()};for(var x in ENV){if(ENV[x]===undefined)delete env[x];else env[x]=ENV[x]}var strings=[];for(var x in env){strings.push(x+"="+env[x])}getEnvStrings.strings=strings}return getEnvStrings.strings}var SYSCALLS={mappings:{},buffers:[null,[],[]],printChar:function(stream,curr){var buffer=SYSCALLS.buffers[stream];if(curr===0||curr===10){(stream===1?out:err)(UTF8ArrayToString(buffer,0));buffer.length=0}else{buffer.push(curr)}},varargs:undefined,get:function(){SYSCALLS.varargs+=4;var ret=HEAP32[SYSCALLS.varargs-4>>2];return ret},getStr:function(ptr){var ret=UTF8ToString(ptr);return ret},get64:function(low,high){return low}};function _environ_get(__environ,environ_buf){var bufSize=0;getEnvStrings().forEach(function(string,i){var ptr=environ_buf+bufSize;HEAP32[__environ+i*4>>2]=ptr;writeAsciiToMemory(string,ptr);bufSize+=string.length+1});return 0}function _environ_sizes_get(penviron_count,penviron_buf_size){var strings=getEnvStrings();HEAP32[penviron_count>>2]=strings.length;var bufSize=0;strings.forEach(function(string){bufSize+=string.length+1});HEAP32[penviron_buf_size>>2]=bufSize;return 0}function _getTempRet0(){return getTempRet0()}function _llvm_eh_typeid_for(type){return type}function __isLeapYear(year){return year%4===0&&(year%100!==0||year%400===0)}function __arraySum(array,index){var sum=0;for(var i=0;i<=index;sum+=array[i++]){}return sum}var __MONTH_DAYS_LEAP=[31,29,31,30,31,30,31,31,30,31,30,31];var __MONTH_DAYS_REGULAR=[31,28,31,30,31,30,31,31,30,31,30,31];function __addDays(date,days){var newDate=new Date(date.getTime());while(days>0){var leap=__isLeapYear(newDate.getFullYear());var currentMonth=newDate.getMonth();var daysInCurrentMonth=(leap?__MONTH_DAYS_LEAP:__MONTH_DAYS_REGULAR)[currentMonth];if(days>daysInCurrentMonth-newDate.getDate()){days-=daysInCurrentMonth-newDate.getDate()+1;newDate.setDate(1);if(currentMonth<11){newDate.setMonth(currentMonth+1)}else{newDate.setMonth(0);newDate.setFullYear(newDate.getFullYear()+1)}}else{newDate.setDate(newDate.getDate()+days);return newDate}}return newDate}function _strftime(s,maxsize,format,tm){var tm_zone=HEAP32[tm+40>>2];var date={tm_sec:HEAP32[tm>>2],tm_min:HEAP32[tm+4>>2],tm_hour:HEAP32[tm+8>>2],tm_mday:HEAP32[tm+12>>2],tm_mon:HEAP32[tm+16>>2],tm_year:HEAP32[tm+20>>2],tm_wday:HEAP32[tm+24>>2],tm_yday:HEAP32[tm+28>>2],tm_isdst:HEAP32[tm+32>>2],tm_gmtoff:HEAP32[tm+36>>2],tm_zone:tm_zone?UTF8ToString(tm_zone):""};var pattern=UTF8ToString(format);var EXPANSION_RULES_1={"%c":"%a %b %d %H:%M:%S %Y","%D":"%m/%d/%y","%F":"%Y-%m-%d","%h":"%b","%r":"%I:%M:%S %p","%R":"%H:%M","%T":"%H:%M:%S","%x":"%m/%d/%y","%X":"%H:%M:%S","%Ec":"%c","%EC":"%C","%Ex":"%m/%d/%y","%EX":"%H:%M:%S","%Ey":"%y","%EY":"%Y","%Od":"%d","%Oe":"%e","%OH":"%H","%OI":"%I","%Om":"%m","%OM":"%M","%OS":"%S","%Ou":"%u","%OU":"%U","%OV":"%V","%Ow":"%w","%OW":"%W","%Oy":"%y"};for(var rule in EXPANSION_RULES_1){pattern=pattern.replace(new RegExp(rule,"g"),EXPANSION_RULES_1[rule])}var WEEKDAYS=["Sunday","Monday","Tuesday","Wednesday","Thursday","Friday","Saturday"];var MONTHS=["January","February","March","April","May","June","July","August","September","October","November","December"];function leadingSomething(value,digits,character){var str=typeof value==="number"?value.toString():value||"";while(str.length<digits){str=character[0]+str}return str}function leadingNulls(value,digits){return leadingSomething(value,digits,"0")}function compareByDay(date1,date2){function sgn(value){return value<0?-1:value>0?1:0}var compare;if((compare=sgn(date1.getFullYear()-date2.getFullYear()))===0){if((compare=sgn(date1.getMonth()-date2.getMonth()))===0){compare=sgn(date1.getDate()-date2.getDate())}}return compare}function getFirstWeekStartDate(janFourth){switch(janFourth.getDay()){case 0:return new Date(janFourth.getFullYear()-1,11,29);case 1:return janFourth;case 2:return new Date(janFourth.getFullYear(),0,3);case 3:return new Date(janFourth.getFullYear(),0,2);case 4:return new Date(janFourth.getFullYear(),0,1);case 5:return new Date(janFourth.getFullYear()-1,11,31);case 6:return new Date(janFourth.getFullYear()-1,11,30)}}function getWeekBasedYear(date){var thisDate=__addDays(new Date(date.tm_year+1900,0,1),date.tm_yday);var janFourthThisYear=new Date(thisDate.getFullYear(),0,4);var janFourthNextYear=new Date(thisDate.getFullYear()+1,0,4);var firstWeekStartThisYear=getFirstWeekStartDate(janFourthThisYear);var firstWeekStartNextYear=getFirstWeekStartDate(janFourthNextYear);if(compareByDay(firstWeekStartThisYear,thisDate)<=0){if(compareByDay(firstWeekStartNextYear,thisDate)<=0){return thisDate.getFullYear()+1}else{return thisDate.getFullYear()}}else{return thisDate.getFullYear()-1}}var EXPANSION_RULES_2={"%a":function(date){return WEEKDAYS[date.tm_wday].substring(0,3)},"%A":function(date){return WEEKDAYS[date.tm_wday]},"%b":function(date){return MONTHS[date.tm_mon].substring(0,3)},"%B":function(date){return MONTHS[date.tm_mon]},"%C":function(date){var year=date.tm_year+1900;return leadingNulls(year/100|0,2)},"%d":function(date){return leadingNulls(date.tm_mday,2)},"%e":function(date){return leadingSomething(date.tm_mday,2," ")},"%g":function(date){return getWeekBasedYear(date).toString().substring(2)},"%G":function(date){return getWeekBasedYear(date)},"%H":function(date){return leadingNulls(date.tm_hour,2)},"%I":function(date){var twelveHour=date.tm_hour;if(twelveHour==0)twelveHour=12;else if(twelveHour>12)twelveHour-=12;return leadingNulls(twelveHour,2)},"%j":function(date){return leadingNulls(date.tm_mday+__arraySum(__isLeapYear(date.tm_year+1900)?__MONTH_DAYS_LEAP:__MONTH_DAYS_REGULAR,date.tm_mon-1),3)},"%m":function(date){return leadingNulls(date.tm_mon+1,2)},"%M":function(date){return leadingNulls(date.tm_min,2)},"%n":function(){return"\n"},"%p":function(date){if(date.tm_hour>=0&&date.tm_hour<12){return"AM"}else{return"PM"}},"%S":function(date){return leadingNulls(date.tm_sec,2)},"%t":function(){return"\t"},"%u":function(date){return date.tm_wday||7},"%U":function(date){var janFirst=new Date(date.tm_year+1900,0,1);var firstSunday=janFirst.getDay()===0?janFirst:__addDays(janFirst,7-janFirst.getDay());var endDate=new Date(date.tm_year+1900,date.tm_mon,date.tm_mday);if(compareByDay(firstSunday,endDate)<0){var februaryFirstUntilEndMonth=__arraySum(__isLeapYear(endDate.getFullYear())?__MONTH_DAYS_LEAP:__MONTH_DAYS_REGULAR,endDate.getMonth()-1)-31;var firstSundayUntilEndJanuary=31-firstSunday.getDate();var days=firstSundayUntilEndJanuary+februaryFirstUntilEndMonth+endDate.getDate();return leadingNulls(Math.ceil(days/7),2)}return compareByDay(firstSunday,janFirst)===0?"01":"00"},"%V":function(date){var janFourthThisYear=new Date(date.tm_year+1900,0,4);var janFourthNextYear=new Date(date.tm_year+1901,0,4);var firstWeekStartThisYear=getFirstWeekStartDate(janFourthThisYear);var firstWeekStartNextYear=getFirstWeekStartDate(janFourthNextYear);var endDate=__addDays(new Date(date.tm_year+1900,0,1),date.tm_yday);if(compareByDay(endDate,firstWeekStartThisYear)<0){return"53"}if(compareByDay(firstWeekStartNextYear,endDate)<=0){return"01"}var daysDifference;if(firstWeekStartThisYear.getFullYear()<date.tm_year+1900){daysDifference=date.tm_yday+32-firstWeekStartThisYear.getDate()}else{daysDifference=date.tm_yday+1-firstWeekStartThisYear.getDate()}return leadingNulls(Math.ceil(daysDifference/7),2)},"%w":function(date){return date.tm_wday},"%W":function(date){var janFirst=new Date(date.tm_year,0,1);var firstMonday=janFirst.getDay()===1?janFirst:__addDays(janFirst,janFirst.getDay()===0?1:7-janFirst.getDay()+1);var endDate=new Date(date.tm_year+1900,date.tm_mon,date.tm_mday);if(compareByDay(firstMonday,endDate)<0){var februaryFirstUntilEndMonth=__arraySum(__isLeapYear(endDate.getFullYear())?__MONTH_DAYS_LEAP:__MONTH_DAYS_REGULAR,endDate.getMonth()-1)-31;var firstMondayUntilEndJanuary=31-firstMonday.getDate();var days=firstMondayUntilEndJanuary+februaryFirstUntilEndMonth+endDate.getDate();return leadingNulls(Math.ceil(days/7),2)}return compareByDay(firstMonday,janFirst)===0?"01":"00"},"%y":function(date){return(date.tm_year+1900).toString().substring(2)},"%Y":function(date){return date.tm_year+1900},"%z":function(date){var off=date.tm_gmtoff;var ahead=off>=0;off=Math.abs(off)/60;off=off/60*100+off%60;return(ahead?"+":"-")+String("0000"+off).slice(-4)},"%Z":function(date){return date.tm_zone},"%%":function(){return"%"}};for(var rule in EXPANSION_RULES_2){if(pattern.includes(rule)){pattern=pattern.replace(new RegExp(rule,"g"),EXPANSION_RULES_2[rule](date))}}var bytes=intArrayFromString(pattern,false);if(bytes.length>maxsize){return 0}writeArrayToMemory(bytes,s);return bytes.length-1}function _strftime_l(s,maxsize,format,tm){return _strftime(s,maxsize,format,tm)}function intArrayFromString(stringy,dontAddNull,length){var len=length>0?length:lengthBytesUTF8(stringy)+1;var u8array=new Array(len);var numBytesWritten=stringToUTF8Array(stringy,u8array,0,u8array.length);if(dontAddNull)u8array.length=numBytesWritten;return u8array}var asmLibraryArg={"x":___assert_fail,"o":___cxa_allocate_exception,"r":___cxa_begin_catch,"s":___cxa_end_catch,"b":___cxa_find_matching_catch_2,"l":___cxa_find_matching_catch_3,"A":___cxa_find_matching_catch_4,"q":___cxa_free_exception,"K":___cxa_rethrow,"p":___cxa_throw,"O":___cxa_uncaught_exceptions,"f":___resumeException,"J":_abort,"R":_emscripten_memcpy_big,"I":_emscripten_resize_heap,"T":_environ_get,"N":_environ_sizes_get,"a":_getTempRet0,"t":invoke_dii,"P":invoke_i,"g":invoke_ii,"B":invoke_iid,"c":invoke_iii,"h":invoke_iiii,"n":invoke_iiiii,"Q":invoke_iiiiid,"C":invoke_iiiiii,"y":invoke_iiiiiii,"L":invoke_iiiiiiii,"E":invoke_iiiiiiiiiiii,"k":invoke_v,"i":invoke_vi,"m":invoke_vid,"G":invoke_vidddddd,"v":invoke_vif,"d":invoke_vii,"e":invoke_viii,"F":invoke_viiidd,"j":invoke_viiii,"H":invoke_viiiii,"M":invoke_viiiiii,"u":invoke_viiiiiii,"z":invoke_viiiiiiiiii,"D":invoke_viiiiiiiiiiiiiii,"w":_llvm_eh_typeid_for,"S":_strftime_l};var asm=createWasm();var ___wasm_call_ctors=Module["___wasm_call_ctors"]=function(){return(___wasm_call_ctors=Module["___wasm_call_ctors"]=Module["asm"]["V"]).apply(null,arguments)};var _pdb2cif=Module["_pdb2cif"]=function(){return(_pdb2cif=Module["_pdb2cif"]=Module["asm"]["W"]).apply(null,arguments)};var _free=Module["_free"]=function(){return(_free=Module["_free"]=Module["asm"]["X"]).apply(null,arguments)};var _mtz2cif=Module["_mtz2cif"]=function(){return(_mtz2cif=Module["_mtz2cif"]=Module["asm"]["Y"]).apply(null,arguments)};var _mxdepo=Module["_mxdepo"]=function(){return(_mxdepo=Module["_mxdepo"]=Module["asm"]["Z"]).apply(null,arguments)};var _cif2pdb=Module["_cif2pdb"]=function(){return(_cif2pdb=Module["_cif2pdb"]=Module["asm"]["_"]).apply(null,arguments)};var _cif2mtz=Module["_cif2mtz"]=function(){return(_cif2mtz=Module["_cif2mtz"]=Module["asm"]["$"]).apply(null,arguments)};var _get_version=Module["_get_version"]=function(){return(_get_version=Module["_get_version"]=Module["asm"]["aa"]).apply(null,arguments)};var _get_str2=Module["_get_str2"]=function(){return(_get_str2=Module["_get_str2"]=Module["asm"]["ba"]).apply(null,arguments)};var _clear_string=Module["_clear_string"]=function(){return(_clear_string=Module["_clear_string"]=Module["asm"]["ca"]).apply(null,arguments)};var _get_global_str_size=Module["_get_global_str_size"]=function(){return(_get_global_str_size=Module["_get_global_str_size"]=Module["asm"]["da"]).apply(null,arguments)};var _malloc=Module["_malloc"]=function(){return(_malloc=Module["_malloc"]=Module["asm"]["ea"]).apply(null,arguments)};var stackSave=Module["stackSave"]=function(){return(stackSave=Module["stackSave"]=Module["asm"]["fa"]).apply(null,arguments)};var stackRestore=Module["stackRestore"]=function(){return(stackRestore=Module["stackRestore"]=Module["asm"]["ga"]).apply(null,arguments)};var stackAlloc=Module["stackAlloc"]=function(){return(stackAlloc=Module["stackAlloc"]=Module["asm"]["ha"]).apply(null,arguments)};var _setThrew=Module["_setThrew"]=function(){return(_setThrew=Module["_setThrew"]=Module["asm"]["ia"]).apply(null,arguments)};var ___cxa_can_catch=Module["___cxa_can_catch"]=function(){return(___cxa_can_catch=Module["___cxa_can_catch"]=Module["asm"]["ja"]).apply(null,arguments)};var ___cxa_is_pointer_type=Module["___cxa_is_pointer_type"]=function(){return(___cxa_is_pointer_type=Module["___cxa_is_pointer_type"]=Module["asm"]["ka"]).apply(null,arguments)};function invoke_ii(index,a1){var sp=stackSave();try{return wasmTable.get(index)(a1)}catch(e){stackRestore(sp);if(e!==e+0&&e!=="longjmp")throw e;_setThrew(1,0)}}function invoke_iii(index,a1,a2){var sp=stackSave();try{return wasmTable.get(index)(a1,a2)}catch(e){stackRestore(sp);if(e!==e+0&&e!=="longjmp")throw e;_setThrew(1,0)}}function invoke_viii(index,a1,a2,a3){var sp=stackSave();try{wasmTable.get(index)(a1,a2,a3)}catch(e){stackRestore(sp);if(e!==e+0&&e!=="longjmp")throw e;_setThrew(1,0)}}function invoke_vii(index,a1,a2){var sp=stackSave();try{wasmTable.get(index)(a1,a2)}catch(e){stackRestore(sp);if(e!==e+0&&e!=="longjmp")throw e;_setThrew(1,0)}}function invoke_vi(index,a1){var sp=stackSave();try{wasmTable.get(index)(a1)}catch(e){stackRestore(sp);if(e!==e+0&&e!=="longjmp")throw e;_setThrew(1,0)}}function invoke_v(index){var sp=stackSave();try{wasmTable.get(index)()}catch(e){stackRestore(sp);if(e!==e+0&&e!=="longjmp")throw e;_setThrew(1,0)}}function invoke_viiii(index,a1,a2,a3,a4){var sp=stackSave();try{wasmTable.get(index)(a1,a2,a3,a4)}catch(e){stackRestore(sp);if(e!==e+0&&e!=="longjmp")throw e;_setThrew(1,0)}}function invoke_iiiii(index,a1,a2,a3,a4){var sp=stackSave();try{return wasmTable.get(index)(a1,a2,a3,a4)}catch(e){stackRestore(sp);if(e!==e+0&&e!=="longjmp")throw e;_setThrew(1,0)}}function invoke_iiii(index,a1,a2,a3){var sp=stackSave();try{return wasmTable.get(index)(a1,a2,a3)}catch(e){stackRestore(sp);if(e!==e+0&&e!=="longjmp")throw e;_setThrew(1,0)}}function invoke_viiiii(index,a1,a2,a3,a4,a5){var sp=stackSave();try{wasmTable.get(index)(a1,a2,a3,a4,a5)}catch(e){stackRestore(sp);if(e!==e+0&&e!=="longjmp")throw e;_setThrew(1,0)}}function invoke_iiiiiii(index,a1,a2,a3,a4,a5,a6){var sp=stackSave();try{return wasmTable.get(index)(a1,a2,a3,a4,a5,a6)}catch(e){stackRestore(sp);if(e!==e+0&&e!=="longjmp")throw e;_setThrew(1,0)}}function invoke_viiiiii(index,a1,a2,a3,a4,a5,a6){var sp=stackSave();try{wasmTable.get(index)(a1,a2,a3,a4,a5,a6)}catch(e){stackRestore(sp);if(e!==e+0&&e!=="longjmp")throw e;_setThrew(1,0)}}function invoke_vidddddd(index,a1,a2,a3,a4,a5,a6,a7){var sp=stackSave();try{wasmTable.get(index)(a1,a2,a3,a4,a5,a6,a7)}catch(e){stackRestore(sp);if(e!==e+0&&e!=="longjmp")throw e;_setThrew(1,0)}}function invoke_iiiiii(index,a1,a2,a3,a4,a5){var sp=stackSave();try{return wasmTable.get(index)(a1,a2,a3,a4,a5)}catch(e){stackRestore(sp);if(e!==e+0&&e!=="longjmp")throw e;_setThrew(1,0)}}function invoke_iid(index,a1,a2){var sp=stackSave();try{return wasmTable.get(index)(a1,a2)}catch(e){stackRestore(sp);if(e!==e+0&&e!=="longjmp")throw e;_setThrew(1,0)}}function invoke_dii(index,a1,a2){var sp=stackSave();try{return wasmTable.get(index)(a1,a2)}catch(e){stackRestore(sp);if(e!==e+0&&e!=="longjmp")throw e;_setThrew(1,0)}}function invoke_viiidd(index,a1,a2,a3,a4,a5){var sp=stackSave();try{wasmTable.get(index)(a1,a2,a3,a4,a5)}catch(e){stackRestore(sp);if(e!==e+0&&e!=="longjmp")throw e;_setThrew(1,0)}}function invoke_viiiiiii(index,a1,a2,a3,a4,a5,a6,a7){var sp=stackSave();try{wasmTable.get(index)(a1,a2,a3,a4,a5,a6,a7)}catch(e){stackRestore(sp);if(e!==e+0&&e!=="longjmp")throw e;_setThrew(1,0)}}function invoke_vid(index,a1,a2){var sp=stackSave();try{wasmTable.get(index)(a1,a2)}catch(e){stackRestore(sp);if(e!==e+0&&e!=="longjmp")throw e;_setThrew(1,0)}}function invoke_vif(index,a1,a2){var sp=stackSave();try{wasmTable.get(index)(a1,a2)}catch(e){stackRestore(sp);if(e!==e+0&&e!=="longjmp")throw e;_setThrew(1,0)}}function invoke_iiiiid(index,a1,a2,a3,a4,a5){var sp=stackSave();try{return wasmTable.get(index)(a1,a2,a3,a4,a5)}catch(e){stackRestore(sp);if(e!==e+0&&e!=="longjmp")throw e;_setThrew(1,0)}}function invoke_iiiiiiii(index,a1,a2,a3,a4,a5,a6,a7){var sp=stackSave();try{return wasmTable.get(index)(a1,a2,a3,a4,a5,a6,a7)}catch(e){stackRestore(sp);if(e!==e+0&&e!=="longjmp")throw e;_setThrew(1,0)}}function invoke_iiiiiiiiiiii(index,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11){var sp=stackSave();try{return wasmTable.get(index)(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11)}catch(e){stackRestore(sp);if(e!==e+0&&e!=="longjmp")throw e;_setThrew(1,0)}}function invoke_viiiiiiiiii(index,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10){var sp=stackSave();try{wasmTable.get(index)(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)}catch(e){stackRestore(sp);if(e!==e+0&&e!=="longjmp")throw e;_setThrew(1,0)}}function invoke_viiiiiiiiiiiiiii(index,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15){var sp=stackSave();try{wasmTable.get(index)(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15)}catch(e){stackRestore(sp);if(e!==e+0&&e!=="longjmp")throw e;_setThrew(1,0)}}function invoke_i(index){var sp=stackSave();try{return wasmTable.get(index)()}catch(e){stackRestore(sp);if(e!==e+0&&e!=="longjmp")throw e;_setThrew(1,0)}}Module["ccall"]=ccall;Module["getValue"]=getValue;Module["UTF8ToString"]=UTF8ToString;Module["writeArrayToMemory"]=writeArrayToMemory;var calledRun;function ExitStatus(status){this.name="ExitStatus";this.message="Program terminated with exit("+status+")";this.status=status}dependenciesFulfilled=function runCaller(){if(!calledRun)run();if(!calledRun)dependenciesFulfilled=runCaller};function run(args){args=args||arguments_;if(runDependencies>0){return}preRun();if(runDependencies>0){return}function doRun(){if(calledRun)return;calledRun=true;Module["calledRun"]=true;if(ABORT)return;initRuntime();if(Module["onRuntimeInitialized"])Module["onRuntimeInitialized"]();postRun()}if(Module["setStatus"]){Module["setStatus"]("Running...");setTimeout(function(){setTimeout(function(){Module["setStatus"]("")},1);doRun()},1)}else{doRun()}}Module["run"]=run;if(Module["preInit"]){if(typeof Module["preInit"]=="function")Module["preInit"]=[Module["preInit"]];while(Module["preInit"].length>0){Module["preInit"].pop()()}}run();
