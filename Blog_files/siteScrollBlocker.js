(window.webpackJsonp__wix_thunderbolt_app=window.webpackJsonp__wix_thunderbolt_app||[]).push([[71],{336:function(e,t,n){"use strict";n.r(t),n.d(t,"site",(function(){return f})),n.d(t,"name",(function(){return r.b})),n.d(t,"SiteScrollBlockerSymbol",(function(){return r.a})),n.d(t,"ISiteScrollBlocker",(function(){return i.ISiteScrollBlocker})),n.d(t,"IScrollBlockedListener",(function(){return i.IScrollBlockedListener}));var r=n(230),i=n(357),o=n(0),c=n(246),s=n(7),u=n(28),l=n(68),a=n.n(l),d=Object(s.h)([u.a],(function(e){var t=0,n=new Map,r=[],i=0,s=function(){return{bodyElement:e.document.body}},u=function(t){var r=s().bodyElement;a.a.measure((function(){i=e.scrollY,r.style.setProperty("--blocked-site-scroll-margin-top",Math.max(.5,i)+"px"),r.classList.add("blockSiteScrolling")})),n.forEach((function(e){var n=e.handleBlockedBy;return n&&n(t)}))},l=function(t){var r=s().bodyElement;r.classList.remove("blockSiteScrolling"),r.style.removeProperty("--blocked-site-scroll-margin-top"),e.scrollTo(0,i),n.forEach((function(e){var n=e.handleUnblockedBy;return n&&n(t)}))};return{setSiteScrollingBlocked:function(t,n){var i;if(!Object(c.e)(e))return t?(i=n,void(1===(r=r.includes(i)?r:Object(o.g)(r,[i])).length&&u(i))):function(e){var t=Object(o.e)(r,1)[0];r=r.filter((function(t){return t!==e}));var n=Object(o.e)(r,1)[0];t!==n&&(n?u(e):l(e))}(n)},registerScrollBlockedListener:function(e){var r=++t;return n.set(r,e),r},unRegisterScrollBlockedListener:function(e){n.delete(e)}}})),f=function(e){e(r.a).to(d)}},357:function(e,t){},68:function(e,t,n){var r;!function(t){"use strict";var i=function(){},o=t.requestAnimationFrame||t.webkitRequestAnimationFrame||t.mozRequestAnimationFrame||t.msRequestAnimationFrame||function(e){return setTimeout(e,16)};function c(){this.reads=[],this.writes=[],this.raf=o.bind(t),i("initialized",this)}function s(e){e.scheduled||(e.scheduled=!0,e.raf(u.bind(null,e)),i("flush scheduled"))}function u(e){i("flush");var t,n=e.writes,r=e.reads;try{i("flushing reads",r.length),l(r),i("flushing writes",n.length),l(n)}catch(e){t=e}if(e.scheduled=!1,(r.length||n.length)&&s(e),t){if(i("task errored",t.message),!e.catch)throw t;e.catch(t)}}function l(e){var t;for(i("run tasks");t=e.shift();)t()}function a(e,t){var n=e.indexOf(t);return!!~n&&!!e.splice(n,1)}c.prototype={constructor:c,measure:function(e,t){i("measure");var n=t?e.bind(t):e;return this.reads.push(n),s(this),n},mutate:function(e,t){i("mutate");var n=t?e.bind(t):e;return this.writes.push(n),s(this),n},clear:function(e){return i("clear",e),a(this.reads,e)||a(this.writes,e)},extend:function(e){if(i("extend",e),"object"!=typeof e)throw new Error("expected object");var t=Object.create(this);return function(e,t){for(var n in t)t.hasOwnProperty(n)&&(e[n]=t[n])}(t,e),t.fastdom=this,t.initialize&&t.initialize(),t},catch:null};var d=t.fastdom=t.fastdom||new c;void 0===(r=function(){return d}.call(d,n,d,e))||(e.exports=r)}("undefined"!=typeof window?window:this)}}]);
//# sourceMappingURL=https://static.parastorage.com/services/wix-thunderbolt/dist/siteScrollBlocker.d442ecb7.chunk.min.js.map